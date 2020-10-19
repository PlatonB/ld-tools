__version__ = 'V1.0'

import os, re, sqlite3

def create_src_dict(src_dir_path, src_file_name, meta_lines_quan, intgen_convdb_path):
        '''
        Функция создаёт словарь {'chrN': [pos, 'rsID'], ...}
        по биаллельным вариантам из исходной таблицы.
        '''
        
        #Значительная часть тулкита привязана к функциональности
        #pysam. Этот парсер биоинформатических таблиц поддерживает
        #запросы по геномным координатам, поэтому дальше будем
        #оперировать в основном ими. Извлечём из конвертационной
        #базы имена хромосом и позиции исходных вариантов, но
        #про rsIDs забывать не будем. Последние пригодятся для
        #планируемого далее уточения идентичности исходного
        #варианта с вариантом из 1000G, т.к. для одной позиции
        #может существовать несколько признанных биаллельными
        #вариантов. Накопим всё в словарь таким образом, чтобы
        #хромосомы стали ключами, а соответствующие остальные
        #элементы - значениями. Поскольку большая таблица БД
        #состоит только из биаллельных вариантов, имеющих rsIDs,
        #неподходящие исходные варианты в словаре не окажутся.
        with open(os.path.join(src_dir_path, src_file_name)) as src_file_opened:
                for meta_line_index in range(meta_lines_quan):
                        src_file_opened.readline()
                rs_ids = set()
                for line in src_file_opened:
                        try:
                                rs_ids.add(re.search(r'rs\d+\b', line).group())
                        except AttributeError:
                                continue
                        
        #Если в исходной таблице
        #не выявилось ни одного
        #валидного варианта, функция
        #возвратит пустой словарь.
        if rs_ids == set():
                return {}
        
        #Преобразование набора rsIDs из множества
        #в кортеж нужно для того, чтобы обеспечить
        #построение SQL-запроса по этому набору.
        #Сам запрос позволяет получить координаты
        #и заново идентификаторы вариантов по rsIDs,
        #пользуясь CHROM-POS-ID-таблицей конвертационной
        #базы данных. Наборы позиций и rsIDs программа
        #сгруппирует по хромосомным ключам. Потом,
        #благодаря такой структуре данных, можно
        #будет безболезненно проводить вычисления
        #для групп вариантов различных хромосом.
        rs_ids, data_by_chrs = tuple(rs_ids), {}
        get_var_basic_info = f'SELECT * FROM variants WHERE ID IN {rs_ids}'.replace(',)', ')')
        with sqlite3.connect(intgen_convdb_path) as conn:
                cursor = conn.cursor()
                for chrom, pos, rs_id in cursor.execute(get_var_basic_info):
                        try:
                                data_by_chrs[chrom].append([pos, rs_id])
                        except KeyError:
                                data_by_chrs[chrom] = [[pos, rs_id]]
                cursor.close()
                
        return data_by_chrs
