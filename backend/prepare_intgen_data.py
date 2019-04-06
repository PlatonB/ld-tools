__version__ = 'V3.0'

import urllib.request, re, os, json, gzip, sys, dbm

def download_intgen_file(intgen_ftp_url, intgen_file_name, intgen_dir_path, index):
        '''
        Скачивание данных 1000 Genomes.
        '''
        
        #Формирование ссылки на скачивание файла 1000
        #Genomes из официального FTP этого проекта, абсолютного
        #пути к этому файлу на компьютере пользователя и часть
        #имени файла без расширения(-ий) (далее - название).
        intgen_file_url = os.path.join(intgen_ftp_url, intgen_file_name)
        intgen_file_path = os.path.join(intgen_dir_path, intgen_file_name)
        intgen_filebase_name = '.'.join(intgen_file_name.split('.')[:index])
        
        #Если 1000 Genomes-файл ещё не скачан, то
        #будут производиться попытки это сделать.
        #В случае, если скачивание срывается из-за
        #сбоев серверов IGSR, эта функция
        #будет рекурсивно запускаться ещё
        #и ещё раз до достижения успеха.
        #После скачивания, или если файл
        #уже есть на диске, функция вернёт
        #путь к нему, а также его название.
        if os.path.exists(intgen_file_path) == False:
                print(f'{intgen_file_name} не найден. Загрузка...')
                try:
                        urllib.request.urlretrieve(intgen_file_url, intgen_file_path)
                        return intgen_file_path, intgen_filebase_name
                except urllib.error.URLError:
                        print('Сбой соединения. Повторная попытка:')
                        return download_intgen_file(intgen_ftp_url, intgen_file_name, intgen_dir_path, index)
        else:
                print(f'{intgen_file_name} уже скачан')
                return intgen_file_path, intgen_filebase_name
        
'''
Отделение от строки символа переноса, её преобразование
в байт-строку с последующей архивацией последней.
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
        SNPs в dbm-базы, из которых можно
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
        
        #Если сохранённой панели нет, скачиваем её.
        #Соответствующая функция также вернёт
        #путь к панели и название таковой.
        print('')
        intgen_samptxt_path, intgen_sampbase_name = download_intgen_file(intgen_hg19_url, intgen_samptxt_name, intgen_dir_path, -1)
        
        #Собираем из названия и расширения имя
        #будущей JSON-редакции панели и путь к ней.
        intgen_sampjson_name = intgen_sampbase_name + '.json'
        intgen_sampjson_path = os.path.join(intgen_dir_path, intgen_sampjson_name)
        
        #JSON-панель создастся только, если она ещё не создана.
        if os.path.exists(intgen_sampjson_path) == True:
                print(f'{intgen_sampjson_name} уже создан\n')
        else:
                print(f'{intgen_sampjson_name} не найден. Создание...\n')
                
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
        
        #Списки для накопления путей к
        #файлам 1000 Genomes без расширений.
        #Расширения будут добавляться
        #уже на уровне фронтендов.
        intgen_vcfbase_paths = []
        
        #В рамках этого цикла будут создаваться и пополняться
        #dbm-базы, основанные на данных из 1000 Genomes-VCFs.
        #Попутно будут скачиваться tabix-индексы.
        for intgen_vcfgz_name in intgen_vcfgz_names:

                #Формирование имени табикс-индекса.
                #При обнаружениии отсутствия этого
                #файла вызываемая функция скачает его.
                intgen_vcftbi_name = intgen_vcfgz_name + '.tbi'
                download_intgen_file(intgen_hg38_url, intgen_vcftbi_name, intgen_dir_path, -3)

                #Проверка на наличие архива 1000 Genomes.
                #Если его нет - он скачается, а также будет
                #возвращён путь к архиву без расширения.
                #Путь будет добавлен в список таких путей.
                intgen_vcfgz_path, intgen_vcfbase_name = download_intgen_file(intgen_hg38_url, intgen_vcfgz_name, intgen_dir_path, -2)
                intgen_vcfbase_paths.append(os.path.join(intgen_dir_path, intgen_vcfbase_name))
                
                #Конструирование имени и пути для далее создаваемой базы.
                intgen_vcfdb_name = intgen_vcfbase_name + '.dbm'
                intgen_vcfdb_path = os.path.join(intgen_dir_path, intgen_vcfdb_name)
                
                #Проверка, нет ли у пользователя уже готовой базы.
                #Если есть, то пропускаем процесс создания.
                if os.path.exists(intgen_vcfdb_path) == True:
                        print(f'{intgen_vcfdb_name} уже создан\n')
                        continue
                print(f'{intgen_vcfdb_name} не найден. Создание...\n')
                
                #Открываем архив 1000 Genomes на чтение таким образом, чтобы 
                #находящийся внутри VCF считывался как строковые данные.
                with gzip.open(intgen_vcfgz_path, mode='rt') as intgen_vcfgz_opened:
                        
                        #Создание dbm-базы.
                        with dbm.open(intgen_vcfdb_path, 'c') as intgen_vcfdb_opened:
                                
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
                                        
                                #Добавление шапки в базу.
                                intgen_vcfdb_opened['header_line'] = compress_line(header_line)
                                
                                #Построчное прочтение основной части VCF-таблицы.
                                for line in intgen_vcfgz_opened:
                                        
                                        #Создание списка, начинающегося
                                        #с идентификатора SNP и
                                        #заканчивающегося ячейкой info-столбца.
                                        #info-ячейка далее пригодится
                                        #для определение аллельности SNP.
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
                                                intgen_vcfdb_opened[part_row[0]] = compress_line(line)
                                                
        return intgen_sampjson_path, intgen_vcfbase_paths
