__version__ = 'V2.0'

import os, urllib.request, sqlite3, re, time
from pysam import tabix_index, VariantFile

def prep_intgen_data(intgen_dir_path):
        '''
        Если даны координаты вариантов, то
        для вычисления LD будет достаточно
        VCF-таблиц проекта 1000 Genomes,
        соответствующих tbi-индексов
        и крохотной таблицы сэмплов.
        Чтобы рассчитывать LD, имея rsIDs,
        потребуется база данных для быстрого
        рандомного доступа к таковым.
        На всех этапах будет выполняться
        проверка, не произведено ли
        то или иное действие ранее.
        '''
        
        #Качаем текстовый файл с
        #именами популяций, гендеров
        #и сэмплов (далее - панель).
        #Он находится в том дистрибутиве
        #1000 Genomes, где используются
        #ещё старые координаты (hg19).
        print('\nsamples.txt', end='... ')
        intgen_samptxt_path = os.path.join(intgen_dir_path, 'samples.txt')
        intgen_hg19page_url = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
        if os.path.exists(intgen_samptxt_path) == False:
                intgen_samptxt_url = os.path.join(intgen_hg19page_url,
                                                  'integrated_call_samples_v3.20130502.ALL.panel')
                urllib.request.urlretrieve(intgen_samptxt_url, intgen_samptxt_path)
        print('OK')
        
        #Инициализация SQLite-БД, необходимой
        #для быстрого получения координат
        #вариантов по идентификаторам, а
        #также сэмплов по полам и популяциям.
        print('conversion.db', end='... ')
        intgen_convdb_path = os.path.join(intgen_dir_path,
                                          'conversion.db')
        conn = sqlite3.connect(intgen_convdb_path)
        cursor = conn.cursor()
        print('OK')
        
        #Верификация структуры панели и размещение её данных в
        #отдельную таблицу базы. Индексация полей этой таблицы
        #не потребуется, поскольку строк в ней очень мало.
        print('samples', end='... ')
        with open(intgen_samptxt_path) as intgen_samptxt_opened:
                header_row = intgen_samptxt_opened.readline().rstrip().split('\t')
                if header_row != ['sample', 'pop', 'super_pop', 'gender']:
                        print('''\nsamples.txt, скорее всего, теперь другой структуры.
Пожалуйста, сообщите об этом в Issues.
Работа программы прервана.''')
                        sys.exit()
                cursor.execute(f'CREATE TABLE IF NOT EXISTS samples ({", ".join(header_row)})')
                cursor.execute('SELECT * FROM samples')
                if cursor.fetchone() is None:
                        samp_data = [line.rstrip().split('\t') for line in intgen_samptxt_opened]
                        cursor.executemany('INSERT INTO samples VALUES (?, ?, ?, ?)', samp_data)
                        conn.commit()
        print('OK')
        
        #Генерация текстового файла
        #со ссылками на размещённые
        #в официальном FTP архивы
        #1000 Genomes. В каждом архиве -
        #VCF-таблица с данными по
        #вариантам одной хромосомы:
        #идентификаторы, некоторые
        #характеристики и, самое
        #главное для вычисления LD -
        #огромный набор фазированных
        #генотипов. Координаты там -
        #hg38. Ссылка на митохондриальную
        #ДНК сохраняться не будет.
        print('urls.txt', end='... ')
        intgen_urlstxt_path = os.path.join(intgen_dir_path, 'urls.txt')
        intgen_hg38page_url = os.path.join(intgen_hg19page_url,
                                           'supporting/GRCh38_positions/')
        if os.path.exists(intgen_urlstxt_path) == False:
                with urllib.request.urlopen(intgen_hg38page_url) as response:
                        intgen_vcf_names = re.findall(r'ALL\.chr(?:\d{1,2}|X|Y)_GRCh38\.genotypes\.\S+?\.vcf\.gz(?=\r\n)',
                                                      response.read().decode('UTF-8'))
                with open(intgen_urlstxt_path, 'w') as intgen_urlstxt_opened:
                        for intgen_vcf_name in intgen_vcf_names:
                                intgen_vcf_url = os.path.join(intgen_hg38page_url,
                                                              intgen_vcf_name)
                                intgen_urlstxt_opened.write(intgen_vcf_url + '\n')
        print('OK')
        
        #При выполнении далее идущего большого блока кода
        #будут скачиваться и индексироваться архивы 1000
        #Genomes, а также пополняться SQLite-таблица
        #соответствия rsIDs с позициями вариантов.
        with open(intgen_urlstxt_path) as intgen_urlstxt_opened:
                for line in intgen_urlstxt_opened:
                        intgen_vcf_url = line.rstrip()
                        
                        #Вытаскиваем из имени архива имя хромосомы.
                        chr_name = re.search(r'(?<=chr)(?:\d{1,2}|X|Y)',
                                             os.path.basename(intgen_vcf_url)).group()
                        
                        #Раскомментируйте и отредактируйте
                        #этот код, если вам нужно будет
                        #работать только с частью хромосом.
##                        if chr_name not in ['6', '14']:
##                                continue
                        
                        #Скачивание архива и его индексация
                        #с помощью tabix. В случае плохого
                        #соединения повторная попытка
                        #загрузки будет производиться через
                        #60 секунд после каждого сбоя. Для
                        #простоты имена архива и индекса
                        #будут состоять только из имени
                        #хромосомы и файловых расширений.
                        print(f'\n{chr_name}.vcf.gz', end='... ')
                        intgen_vcf_path = os.path.join(intgen_dir_path,
                                                       f'{chr_name}.vcf.gz')
                        if os.path.exists(intgen_vcf_path) == False:
                                while True:
                                        try:
                                                urllib.request.urlretrieve(intgen_vcf_url,
                                                                           intgen_vcf_path)
                                                break
                                        except:
                                                print('Fail')
                                                os.remove(intgen_vcf_path)
                                                time.sleep(60)
                                                print(f'{chr_name}.vcf.gz', end='... ')
                        print('OK')
                        print(f'{chr_name}.vcf.gz.tbi', end='... ')
                        if os.path.exists(intgen_vcf_path + '.tbi') == False:
                                tabix_index(intgen_vcf_path, preset='vcf')
                        print('OK')
                        
                        #Создание таблицы соответствия.
                        print('variants', end='... ')
                        cursor.execute('CREATE TABLE IF NOT EXISTS variants (CHROM TEXT, POS INTEGER, ID TEXT)')
                        cursor.execute(f'SELECT * FROM variants WHERE CHROM = "{chr_name}"')
                        if cursor.fetchone() is not None:
                                print('OK')
                                continue
                        
                        #Пополнение таблицы соответствия фрагментами по 100000
                        #строк. Варианты с не-rs-идентификаторами, а также
                        #мультиаллельные варианты в эту таблицу не пойдут, т.к.
                        #для них невозможно или крайне затруднительно считать
                        #LD. Особо коварная разновидность мультиаллельных
                        #вариантов - повторы с разным количеством звеньев.
                        #В 1000 Genomes они оформлены как наборы биаллельных
                        #вариантов. В пределах каждого из наборов - один и тот
                        #же rsID. Программа их распознаёт и тоже выбрасывает.
                        with VariantFile(intgen_vcf_path) as intgen_vcf_opened:
                                fragment, rs_ids, fragment_len, max_fragment_len = [], set(), 0, 100000
                                for rec in intgen_vcf_opened.fetch():
                                        if re.match(r'rs\d+$', rec.id) is None \
                                           or 'MULTI_ALLELIC' in rec.info:
                                                continue
                                        if rec.id in rs_ids:
                                                continue
                                        fragment.append([rec.chrom, rec.pos, rec.id])
                                        rs_ids.add(rec.id)
                                        fragment_len += 1
                                        if fragment_len == max_fragment_len:
                                                cursor.executemany('INSERT INTO variants VALUES (?, ?, ?)', fragment)
                                                fragment.clear()
                                                rs_ids.clear()
                                                fragment_len = 0
                        if fragment_len > 0:
                                cursor.executemany('INSERT INTO variants VALUES (?, ?, ?)', fragment)
                        conn.commit()
                        print('OK')
                        
        #Индексация таблицы соответствия по полю с rsIDs.
        print('\nid', end='... ')
        cursor.execute('CREATE INDEX IF NOT EXISTS "id" ON variants (ID)')
        conn.commit()
        print('OK')
        
        #Дисконнект.
        cursor.close()
        conn.close()
        
        return intgen_convdb_path
