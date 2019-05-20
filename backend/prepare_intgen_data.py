__version__ = 'V4.0'

import urllib.request, re, os, json, gzip, sys, plyvel

def download_intgen_file(intgen_ftp_url, intgen_file_name, intgen_dir_path):
        '''
        Скачивание данных 1000 Genomes.
        '''
        
        #Формирование ссылки на скачивание файла 1000
        #Genomes из официального FTP этого проекта и абсолютного
        #пути к этому файлу на компьютере пользователя.
        intgen_file_url = os.path.join(intgen_ftp_url, intgen_file_name)
        intgen_file_path = os.path.join(intgen_dir_path, intgen_file_name)
        
        #Если 1000 Genomes-файл ещё не скачан, то
        #будут производиться попытки это сделать.
        #В случае, если скачивание срывается
        #из-за сбоев серверов IGSR, эта функция
        #станет рекурсивно запускаться ещё
        #и ещё раз до достижения успеха.
        #После скачивания, или если файл уже
        #есть на диске, функция вернёт путь к нему.
        if os.path.exists(intgen_file_path) == False:
                print(f'{intgen_file_name} не найден. Загрузка...')
                try:
                        urllib.request.urlretrieve(intgen_file_url, intgen_file_path)
                        return intgen_file_path
                except urllib.error.URLError:
                        print('Сбой соединения. Повторная попытка:')
                        return download_intgen_file(intgen_ftp_url, intgen_file_name, intgen_dir_path)
        else:
                print(f'{intgen_file_name} уже скачан')
                return intgen_file_path
        
'''
Отделение от строки символа переноса, её
преобразование в байт-строку, архивация последней.
Искав модуль архивации, я остановился на gzip, и стал
применять его со степенью компрессии по умолчанию (9).
Экспериментальным путём определил, что именно такой
выбор оптимален по соотношению скорость-сжатие.
'''
compress_line = lambda line: gzip.compress(line.split('\n')[0].encode())

def process_intgen_data(intgen_dir_path):
        '''
        Формирование JSON-редакции панели сэмплов.
        Размещение данных 1000 Genomes по
        SNPs в LevelDB-базу, из которой можно
        будет быстро извлекать строки
        по идентификаторам SNP.
        '''
        
        #Если папка для 1000 Genomes-данных и
        #соответствующих баз абсолютно пуста, то
        #всё будет скачиваться и создаваться с нуля.
        #Поскольку это - процесс долгий, то пользователю
        #будет показано соответствующее предупреждение.
        if os.listdir(intgen_dir_path) == []:
                print('''\nДля рассчёта LD необходимы данные 1000 Genomes.
Их загрузка и оптимизация могут занять сутки.''')
                
        #Захардкодим URLs тех каталогов FTP IGSR,
        #где хостятся актуальные данные, необходимые
        #для рассчёта LD (1000 Genomes Phase 3).
        #В одном каталоге координаты соответствуют
        #сборке GRCh37 (hg19), в другом - GRCh38 (hg38).
        #Из первого нужно будет скачать текстовый файл
        #с популяциями и сэмплами (далее - панель).
        #Из второго - архивы 1000 Genomes и их индексы.
        #В каждом архиве - VCF-таблица с данными по SNPs
        #одной хромосомы: идентификатор снипа, некоторые
        #его характеристики и набор фазированных генотипов.
        intgen_hg19_url = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
        intgen_hg38_url = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/'
        
        #Преобразуем страницы обоих каталогов в строки.
        #Из этих строк можно будет извлекать
        #регулярками имена нужных файлов.
        with urllib.request.urlopen(intgen_hg19_url) as response:
                intgen_hg19_page = str(response.read())
        with urllib.request.urlopen(intgen_hg38_url) as response:
                intgen_hg38_page = str(response.read())
                
        #Получаем имя панели.
        intgen_samptxt_name = re.search(r'integrated_call_samples\S+?\.ALL\.panel(?=\\)', intgen_hg19_page).group()
        
        print('')
        
        #Если сохранённой панели нет, скачиваем её.
        #Соответствующая функция также вернёт путь к панели.
        intgen_samptxt_path = download_intgen_file(intgen_hg19_url, intgen_samptxt_name, intgen_dir_path)
        
        #Собираем имя будущей JSON-редакции панели и путь к ней.
        intgen_sampjson_name = '.'.join(intgen_samptxt_name.split('.')[:-1]) + '.json'
        intgen_sampjson_path = os.path.join(intgen_dir_path, intgen_sampjson_name)
        
        #JSON-панель создастся только, если она ещё не создана.
        if os.path.exists(intgen_sampjson_path) == True:
                print(f'{intgen_sampjson_name} уже создан')
        else:
                print(f'{intgen_sampjson_name} не найден. Создание...')
                
                #Открытие исходной (txt) версии панели на чтение.
                with open(intgen_samptxt_path) as intgen_samptxt_opened:
                        
                        #Создание JSON-панели.
                        with open(intgen_sampjson_path, 'w') as intgen_sampjson_opened:
                                
                                #Считывание хэдера, его преобразование в
                                #список и проверка на наличие и соблюдение
                                #последовательности определённых элементов.
                                header_row = intgen_samptxt_opened.readline().split('\n')[0].split('\t')
                                if header_row[:4] != ['sample',
                                                      'pop',
                                                      'super_pop',
                                                      'gender']:
                                        print(f'''{intgen_samptxt_name}, скорее всего, теперь другой структуры.
Пожалуйста, сообщите об этом в Issues.
Работа программы прервана.''')
                                        sys.exit()
                                        
                                #Создание объекта, состоящего из словарей,
                                #вложенных друг в друга в таком порядке:
                                #{суперпопуляция: {субпопуляция: {пол: сэмпл}}}.
                                intgen_samptxt_dict = {}
                                for line in intgen_samptxt_opened:
                                        row = line.split('\n')[0].split('\t')
                                        if row[2] not in intgen_samptxt_dict:
                                                intgen_samptxt_dict[row[2]] = {row[1]: {row[3]: [row[0]]}}
                                        elif row[1] not in intgen_samptxt_dict[row[2]]:
                                                intgen_samptxt_dict[row[2]][row[1]] = {row[3]: [row[0]]}
                                        elif row[3] not in intgen_samptxt_dict[row[2]][row[1]]:
                                                intgen_samptxt_dict[row[2]][row[1]][row[3]] = [row[0]]
                                        else:
                                                intgen_samptxt_dict[row[2]][row[1]][row[3]].append(row[0])
                                                
                                #Прописывание полученного объекта в файл как JSON.
                                json.dump(intgen_samptxt_dict, intgen_sampjson_opened, indent=4)
                                
        #Вытаскиваем имена архивов 1000 Genomes.
        #Данные по митохондриальной ДНК игнорируются.
        intgen_vcfgz_names = re.findall(r'ALL\.chr(?:\d{1,2}|X|Y)_GRCh38\.genotypes\.\S+?\.vcf\.gz(?=\\)', intgen_hg38_page)
        
        #Список для накопления путей к
        #vcf.gz-архивам 1000 Genomes.
        intgen_vcfgz_paths = []
        
        #Создание имени и пути для LevelDB-базы.
        intgen_vcfdb_name = '1000_Genomes_ldb'
        intgen_vcfdb_path = os.path.join(intgen_dir_path, intgen_vcfdb_name)

        print('')
        
        #Проверка, нет ли у пользователя уже готовой БД.
        #Если есть, то пропускаем процесс создания.
        if os.path.exists(intgen_vcfdb_path) == True:
                intgen_vcfdb_exists = True
                print(f'{intgen_vcfdb_name} уже создана')
        else:
                intgen_vcfdb_exists = False
                print(f'{intgen_vcfdb_name} не найдена. Будет создана')
                intgen_vcfdb_opened = plyvel.DB(intgen_vcfdb_path, create_if_missing=True)
                
        #В рамках этого цикла будут скачиваться и/или
        #создаваться все или только недостающие файлы.
        #Это - 1000 Genomes-архивы и tbi-индексы к ним,
        #а также LevelDB-база, которая должна в конечном
        #итоге включать в себя все данные 1000 Genomes.
        for intgen_vcfgz_name in intgen_vcfgz_names:
                
                print('')
                
                #Формирование имени табикс-индекса.
                #При обнаружениии отсутствия этого
                #файла вызываемая функция скачает его.
                intgen_vcftbi_name = intgen_vcfgz_name + '.tbi'
                download_intgen_file(intgen_hg38_url, intgen_vcftbi_name, intgen_dir_path)
                
                #Проверка на наличие архива 1000 Genomes.
                #Если его нет - он скачается.
                #В любом случае, будет возвращён путь к нему.
                #Путь будет добавлен в список таких путей.
                intgen_vcfgz_path = download_intgen_file(intgen_hg38_url, intgen_vcfgz_name, intgen_dir_path)
                intgen_vcfgz_paths.append(intgen_vcfgz_path)

                #Если БД уже создана, то в данном
                #цикле будет выполняться работа только
                #с самими 1000 Genomes-архивами.
                if intgen_vcfdb_exists == True:
                        continue
                
                print(f'{intgen_vcfdb_name} пополняется...')

                #Дальнейшие действия - для
                #пополнения LevelDB-базы.
                #Открываем архив 1000 Genomes
                #на чтение таким образом, чтобы
                #находящийся внутри него VCF
                #считывался как строковые данные.
                with gzip.open(intgen_vcfgz_path, mode='rt') as intgen_vcfgz_opened:
                        
                        #Прочтение вхолостую строк метаинформации.
                        #Получение шапки VCF-таблицы.
                        meta_line = intgen_vcfgz_opened.readline()
                        while meta_line.startswith('##'):
                                meta_line = intgen_vcfgz_opened.readline()
                        header_line = meta_line
                        
                        #Валидация первых 9 элементов шапки.
                        if header_line.split('\t')[:9] != ['#CHROM',
                                                           'POS',
                                                           'ID',
                                                           'REF',
                                                           'ALT',
                                                           'QUAL',
                                                           'FILTER',
                                                           'INFO',
                                                           'FORMAT']:
                                print(f'''{intgen_vcfgz_name}, скорее всего, теперь другой структуры.
Пожалуйста, сообщите об этом в Issues.
Работа программы прервана.''')
                                sys.exit()
                                
                        #Если шапка не добавлена в базу - добавляем.
                        if intgen_vcfdb_opened.get('header_line'.encode()) == None:
                                intgen_vcfdb_opened.put('header_line'.encode(), compress_line(header_line))
                                
                        #Построчное прочтение основной части VCF-таблицы.
                        for line in intgen_vcfgz_opened:
                                
                                #Создание списка, начинающегося
                                #с идентификатора SNP и
                                #заканчивающегося ячейкой info-столбца.
                                #info-ячейка далее пригодится
                                #для определения аллельности SNP.
                                part_row = line.split('\t')[2:8]
                                
                                #Вместо refSNPID могут попадаться
                                #другие типы идентификаторов.
                                #Текущая версия LD-бэкенда
                                #отправляет в базы только те
                                #строки, которые включают refSNPID.
                                #refSNPID становятся ключами,
                                #а содержащие их строки (в
                                #сжатом виде) - значениями.
                                #Мультиаллельные SNP пока не
                                #поддерживаются калькулятором LD,
                                #поэтому соответствующие refSNPID
                                #в базу добавляться не будут.
                                if re.match(r'rs\d+$', part_row[0]) != None and part_row[5].find('MULTI_ALLELIC') == -1:
                                        intgen_vcfdb_opened.put(part_row[0].encode(), compress_line(line))
                                        
        if intgen_vcfdb_exists == False:
                intgen_vcfdb_opened.close()
                
        return intgen_sampjson_path, intgen_vcfdb_path, intgen_vcfgz_paths
