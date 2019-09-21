__version__ = 'V5.0'

import urllib.request, re, os, json, gzip, sys, plyvel

def download_intgen_file(intgen_page_url, intgen_file_name, intgen_dir_path):
        '''
        Скачивание данных 1000 Genomes,
        если они не были скачаны до этого.
        Необходимо для дальнейшего функционирования
        программы в отсутствие интернета.
        '''
        
        #Формирование ссылки на скачивание файла 1000 Genomes
        #из официального FTP этого проекта и абсолютного
        #пути к данному файлу на компьютере пользователя.
        intgen_file_url = os.path.join(intgen_page_url, intgen_file_name)
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
                        return download_intgen_file(intgen_page_url, intgen_file_name, intgen_dir_path)
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
        SNPs в LevelDB-базу, из которой
        можно будет быстро извлекать
        строки по идентификаторам SNP.
        '''
        
        #Если папка для 1000 Genomes-архивов
        #и базы данных абсолютно пуста, то всё
        #будет скачиваться и создаваться с нуля.
        #Поскольку это - процесс долгий, то пользователю
        #будет показано соответствующее предупреждение.
        if os.listdir(intgen_dir_path) == []:
                print('''\nДля рассчёта LD необходимы данные 1000 Genomes.
Их загрузка и оптимизация могут занять сутки.\n''')
                
        #Создаём переменную для
        #ссылки на FTP-каталог,
        #содержащий текстовый
        #файл с популяциями и
        #сэмплами (далее - панель).
        #Другая переменная
        #будет для имени панели.
        intgen_hg19page_url = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
        intgen_samptxt_name = 'integrated_call_samples_v3.20130502.ALL.panel'

        print('')
        
        #Если сохранённой панели нет, скачиваем её.
        #Эта функция также вернёт путь к панели.
        intgen_samptxt_path = download_intgen_file(intgen_hg19page_url, intgen_samptxt_name, intgen_dir_path)
        
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
                                
        #Захардкодим URL каталога FTP IGSR,
        #содержащего данные 1000 Genomes,
        #характеризующиеся обновлёнными
        #до GRCh38 (hg38) координатами.
        #Придумаем имя файла, предназначенного
        #для хранения текста, содержащего
        #имена файлов этого каталога.
        #Потом, если бэкенд не запускался
        #ранее, программа создаст упомянутый
        #файл (далее - hg38-файл) и заполнит
        #его текстовой копией FTP-каталога.
        #Из hg38-файла далее потребуется брать
        #имена только включающих в себя генотипы
        #архивов и соответствующих tabix-индексов.
        #В каждом таком архиве - VCF-таблица
        #с данными по SNPs одной хромосомы:
        #идентификаторы, некоторые сниповые
        #характеристики и, самое главное
        #для вычисления LD - огромный
        #набор фазированных генотипов.
        intgen_hg38page_url = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/'
        intgen_hg38file_name = 'hg38_filenames.txt'
        
        #Сохранение FTP-каталога с именами
        #файлов 1000 Genomes сборки hg38,
        #если это ранее не было сделано.
        #Сохранить его важно для
        #обеспечения последующих запусков
        #программы полностью оффлайн.
        intgen_hg38file_path = os.path.join(intgen_dir_path, intgen_hg38file_name)
        if os.path.exists(intgen_hg38file_path) == True:
                print(f'\n{intgen_hg38file_name} уже создан')
        else:
                print(f'\n{intgen_hg38file_name} не найден. Создание...')
                with urllib.request.urlopen(intgen_hg38page_url) as response:
                        with open(intgen_hg38file_path, 'w') as intgen_hg38file_opened:
                                intgen_hg38file_opened.write(response.read().decode('UTF-8'))
                                
        #Извлекаем текст, содержащий
        #имена hg38-1000Genomes-файлов.
        #Теперь всё готово для отбора
        #регуляркой нужных из них.
        with open(intgen_hg38file_path) as intgen_hg38file_opened:
                hg38_filenames_raw_str = intgen_hg38file_opened.read()
                
        #Вытаскиваем имена архивов 1000 Genomes,
        #несущих данные, обязательные для вычисления LD.
        #Всё, что связано с митохондриальной ДНК, игнорируется.
        intgen_vcfgz_names = re.findall(r'ALL\.chr(?:\d{1,2}|X|Y)_GRCh38\.genotypes\.\S+?\.vcf\.gz(?=\n)', hg38_filenames_raw_str)
        
        #Список для накопления путей к подлежащим
        #скачиванию vcf.gz-архивам 1000 Genomes.
        intgen_vcfgz_paths = []
        
        #Создание имени и пути для LevelDB-базы.
        intgen_vcfdb_name = '1000_Genomes_ldb'
        intgen_vcfdb_path = os.path.join(intgen_dir_path, intgen_vcfdb_name)
        
        #Проверка, нет ли у пользователя уже готовой БД.
        #Если есть, то пропускаем процесс создания.
        if os.path.exists(intgen_vcfdb_path) == True:
                intgen_vcfdb_exists = True
                print(f'\n{intgen_vcfdb_name} уже создана')
        else:
                intgen_vcfdb_exists = False
                print(f'\n{intgen_vcfdb_name} не найдена. Будет создана')
                intgen_vcfdb_opened = plyvel.DB(intgen_vcfdb_path, create_if_missing=True)
                
        #В рамках этого цикла будут
        #скачиваться и/или создаваться
        #все или только недостающие файлы.
        #Это - подробно описанные выше
        #1000 Genomes-архивы и tbi-индексы
        #к ним, а также LevelDB-база,
        #которая должна в конечном итоге
        #включать в себя данные 1000
        #Genomes, позволяющие считать LD.
        for intgen_vcfgz_name in intgen_vcfgz_names:
                
                print('')
                
                #Формирование имени табикс-индекса.
                #При обнаружениии отсутствия этого
                #файла вызываемая функция скачает его.
                intgen_vcftbi_name = intgen_vcfgz_name + '.tbi'
                download_intgen_file(intgen_hg38page_url, intgen_vcftbi_name, intgen_dir_path)
                
                #Проверка на наличие архива 1000 Genomes.
                #Если его нет - он скачается.
                #В любом случае, будет возвращён путь к нему.
                #Путь будет добавлен в список таких путей.
                intgen_vcfgz_path = download_intgen_file(intgen_hg38page_url, intgen_vcfgz_name, intgen_dir_path)
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
