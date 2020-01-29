__version__ = 'V1.0'

import sys, os, urllib.request, re, gzip
sys.dont_write_bytecode = True
from backend.init_intgen_db import *
from pymongo import ASCENDING

def create_intgen_db():
        '''
        Размещение данных 1000 Genomes в MongoDB-базу.
        Из неё можно будет быстро извлекать фазированные
        генотипы по идентификаторам или позициям SNPs.
        В базу также будет включена коллекция сэмплов,
        позволяющая фильтровать фазированные генотипы.
        '''
        
        print('''\nДля рассчёта LD необходимы данные 1000 Genomes.
Их загрузка и размещение в БД могут занять несколько суток.
БД будет весить около 450 ГБ.\n''')
        
        #Создание папки, в которую будут
        #скачиваться файлы 1000 Genomes.
        #Работа с этой папкой далее будет
        #устроена так: каждый скачанный файл
        #просуществует лишь то время, пока
        #данные из него заливаются в базу.
        intgen_tempdir_path = os.path.join(os.getcwd(), '1000_genomes_temp')
        os.mkdir(intgen_tempdir_path)
        
        #Текстовый файл с именами популяций,
        #гендеров и сэмплов (далее - панель)
        #находится в том дистрибутиве
        #1000 Genomes, где используются
        #ещё старые координаты (hg19).
        #Находим там этот файл и качаем.
        intgen_hg19page_url = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
        intgen_samptxt_name = 'integrated_call_samples_v3.20130502.ALL.panel'
        intgen_samptxt_url = os.path.join(intgen_hg19page_url, intgen_samptxt_name)
        intgen_samptxt_path = os.path.join(intgen_tempdir_path, intgen_samptxt_name)
        urllib.request.urlretrieve(intgen_samptxt_url, intgen_samptxt_path)
        
        #Открытие панели на чтение.
        with open(intgen_samptxt_path) as intgen_samptxt_opened:
                
                #Считывание хэдера, его преобразование в
                #список и проверка на наличие и соблюдение
                #последовательности всех элементов.
                header_row = intgen_samptxt_opened.readline().rstrip().split('\t')
                if header_row != ['sample',
                                  'pop',
                                  'super_pop',
                                  'gender']:
                        print(f'''\n{intgen_samptxt_name}, скорее всего, теперь другой структуры.
Пожалуйста, сообщите об этом в Issues.
Работа программы прервана.''')
                        sys.exit()
                        
                print(f'{intgen_samptxt_name} размещается в БД')
                
                #Имена полей этой коллекции
                #без изменений перейдут из шапки.
                #Индексация полей не потребуется,
                #поскольку объём данных очень мал.
                intgen_sampcoll_obj.insert_many([dict(zip(header_row,
                                                          line.rstrip().split('\t'))) for line in intgen_samptxt_opened])
                
        #TXT-версия панели больше не нужна.
        os.remove(intgen_samptxt_path)
        
        #Далее будем работать с основной
        #частью данных 1000 Genomes.
        #Они программными методами
        #приведены к hg38-координатам.
        #Нас будут интересовать только
        #громоздкие gzip-архивы, без
        #соответствующих tabix-индексов.
        #Последние не нужны, т.к. у MongoDB
        #с запасом хватает функционала
        #быстрого извлечения данных.
        #В каждом архиве - VCF-таблица с
        #данными по SNPs одной хромосомы:
        #идентификаторы, некоторые сниповые
        #характеристики и, самое главное
        #для вычисления LD - огромный
        #набор фазированных генотипов.
        intgen_hg38page_url = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/'
        with urllib.request.urlopen(intgen_hg38page_url) as response:
                hg38_filenames_raw_str = response.read().decode('UTF-8')
                
        #Вытаскиваем имена архивов 1000 Genomes,
        #несущих данные, обязательные для вычисления LD.
        #Всё, что связано с митохондриальной ДНК, игнорируется.
        intgen_vcf_names = re.findall(r'ALL\.chr(?:\d{1,2}|X|Y)_GRCh38\.genotypes\.\S+?\.vcf\.gz(?=\r\n)',
                                      hg38_filenames_raw_str)
        
        #В рамках этого цикла будут
        #скачиваться позволяющие
        #считать LD архивы, а данные
        #каждого архива - идти в
        #соответствующую коллекцию.
        #Коллекции для простоты будут
        #называться chr1, chr2 и т.д..
        for intgen_vcf_name in intgen_vcf_names:
                
                #Вытаскиваем из имени архива имя хромосомы.
                chr_name = re.search(r'chr(?:\d{1,2}|X|Y)',
                                     intgen_vcf_name).group()
                
                #Раскомментируйте и отредактируйте
                #этот код, если нужно разместить
                #в базу только часть хромосом.
                #С неполной БД фронтенды будут
                #работать корректно, но выдавать
                #результаты только по снипам
                #ограниченного набора хромосом.
##                if chr_name not in ['chr6', 'chr14']:
##                        continue
                
                print(f'\n{intgen_vcf_name} скачивается')
                
                #Скачивание архива.
                intgen_vcf_url = os.path.join(intgen_hg38page_url, intgen_vcf_name)
                intgen_vcf_path = os.path.join(intgen_tempdir_path, intgen_vcf_name)
                urllib.request.urlretrieve(intgen_vcf_url, intgen_vcf_path)
                
                #Дальнейшие действия - для
                #пополнения MongoDB-базы.
                #Открываем архив 1000 Genomes
                #на чтение таким образом, чтобы
                #находящийся внутри него VCF
                #считывался как строковые данные.
                #Скипаем строки метаинформации.
                with gzip.open(intgen_vcf_path, mode='rt') as intgen_vcf_opened:
                        for line in intgen_vcf_opened:
                                if line.startswith('##') == False:
                                        header_row = line.rstrip().split('\t')
                                        break
                                
                        #Валидация первых 9 элементов шапки.
                        if header_row[:9] != ['#CHROM',
                                              'POS',
                                              'ID',
                                              'REF',
                                              'ALT',
                                              'QUAL',
                                              'FILTER',
                                              'INFO',
                                              'FORMAT']:
                                print(f'''\n{intgen_vcf_name}, скорее всего, теперь другой структуры.
Пожалуйста, сообщите об этом в Issues.
Работа программы прервана.''')
                                sys.exit()
                                
                        #Создание MongoDB-коллекции, предназначенной
                        #для данных только что скачанного архива.
                        intgen_vcfcoll_obj = intgen_db_obj.create_collection(chr_name,
                                                                             storageEngine={'wiredTiger':
                                                                                            {'configString':
                                                                                             'block_compressor=zstd'}})
                        
                        print(f'{intgen_vcf_name} размещается в БД')
                        
                        #Коллекция будет пополняться
                        #порциями по 10000 строк.
                        #Строки - здоровенные, поэтому
                        #лить ещё большими кусками не стоит.
                        #Чревато переполнением RAM.
                        fragment, fragment_len, max_fragment_len = [], 0, 10000
                        
                        #Считывание строк основной части
                        #сжатого VCF с преобразованием в список.
                        for line in intgen_vcf_opened:
                                row = line.rstrip().split('\t')
                                
                                #Вместо refSNPID могут попадаться
                                #другие типы идентификаторов.
                                #Текущая версия LD-бэкенда
                                #отправляет в базы только те
                                #строки, которые включают refSNPID.
                                #Мультиаллельные SNP пока не
                                #поддерживаются калькулятором LD,
                                #поэтому соответствующие refSNPID
                                #в базу добавляться не будут.
                                if re.match(r'rs\d+$', row[2]) != None and \
                                   row[7].find('MULTI_ALLELIC') == -1:
                                        
                                        #MongoDB позволяет добавлять в
                                        #документы вложенные структуры.
                                        #Воспользовавшись этим, сделаем
                                        #более парсибельным INFO-столбец.
                                        #В обычных VCF он представляет
                                        #собой длинную строку, состоящую
                                        #из парных и одиночных элементов.
                                        #Последние называются флагами.
                                        #INFO-подструктура будет словарём:
                                        #парные элементы станут ключами и
                                        #значениями, а флаги пойдут в список,
                                        #присоединённый к отдельному ключу.
                                        info_two_dim = [element.split('=') for element in row[7].split(';')]
                                        info_dict = {'flags': []}
                                        for pair in info_two_dim:
                                                if len(pair) == 2:
                                                        info_dict[pair[0]] = pair[1]
                                                else:
                                                        info_dict['flags'].append(pair[0])
                                                        
                                        #Будущий MongoDB-документ - словарь, создаваемый
                                        #из всех элементов хэдера и текущей строки VCF.
                                        #При этом хэдер мы никак не модифицируем,
                                        #а в строке производим пару изменений:
                                        #задаём координате SNP числовой тип данных,
                                        #а INFO-строку меняем на INFO-словарь.
                                        row[1], row[7] = int(row[1]), info_dict
                                        fragment.append(dict(zip(header_row, row)))
                                        
                                        #Порциональное пополнение MongoDB-коллекции.
                                        fragment_len += 1
                                        if fragment_len == max_fragment_len:
                                                intgen_vcfcoll_obj.insert_many(fragment)
                                                fragment.clear()
                                                fragment_len = 0
                if fragment_len > 0:
                        intgen_vcfcoll_obj.insert_many(fragment)
                        
                print('POS и ID индексируются')
                
                #Индексация полей с координатами и refSNPIDs.
                intgen_vcfcoll_obj.create_index([(header_row[1], ASCENDING)])
                intgen_vcfcoll_obj.create_index([(header_row[2], ASCENDING)])
                
                #VCF залит в базу и
                #больше не пригодится.
                os.remove(intgen_vcf_path)
                
        #Создание базы завершилось,
        #поэтому папку для данных
        #1000 Genomes можно удалить.
        os.rmdir(intgen_tempdir_path)
