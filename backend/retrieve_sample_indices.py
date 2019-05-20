__version__ = 'V1.2'

import json, plyvel, gzip, sys

def retrieve_sample_indices(intgen_sampjson_path, populations, genders, intgen_vcfdb_path):
        '''
        Получение сэмплов, соответствующих заданным популяциям и полам.
        '''
        
        #Открытие JSON-варианта панели сэмплов на чтение.
        #Извлечение хранящегося там объекта в оперативную память.
        with open(intgen_sampjson_path) as intgen_sampjson_opened:
                intgen_sampdict = json.load(intgen_sampjson_opened)
                
        #Открываем ранее созданную
        #1000 Genomes-БД и достаём оттуда
        #шапку, содержащую названия сэмплов.
        #Разархивируем её, конвертируем
        #в Юникод и преобразуем в список.
        intgen_vcfdb_opened = plyvel.DB(intgen_vcfdb_path, create_if_missing=False)
        header_row = gzip.decompress(intgen_vcfdb_opened.get('header_line'.encode())).decode('utf-8').split('\t')
        intgen_vcfdb_opened.close()
        
        #Пользователь может по ошибке указать популяцию дважды или
        #выбрать суперпопуляцию и принадлежащую ей субпопуляцию вместе.
        #В обоих случаях от дублей спасёт размещение сэмплов во множество.
        #Во втором случае суперпопуляция как бы поглотит субпопуляцию.
        samples = set()
        
        #Добавление сэмплов всех популяций в конечное множество
        #при соответствующем пользовательском выборе.
        #Гендеры при этом могут быть как один, так и оба, в
        #зависимости от указанной пользователем опции.
        if populations == ['ALL']:
                for sup_pop in intgen_sampdict:
                        for sub_pop in intgen_sampdict[sup_pop]:
                                for gender in genders:
                                        for sample in intgen_sampdict[sup_pop][sub_pop][gender]:
                                                samples.add(sample)
                                                
        #Пользователь обозначил не все,
        #а только определённые популяции.
        else:
                for population in populations:
                        
                        #Поиск выбранной популяции среди суперпопуляций,
                        #добавление соответствующих сэмплов во множество.
                        if population in intgen_sampdict:
                                for sub_pop in intgen_sampdict[population]:
                                        for gender in genders:
                                                for sample in intgen_sampdict[population][sub_pop][gender]:
                                                        samples.add(sample)
                                                        
                        #Не нашли среди суперпопуляций, тогда ищем
                        #среди суб- и пополняем сэмплами множество.
                        else:
                                for sup_pop in intgen_sampdict:
                                        if population in intgen_sampdict[sup_pop]:
                                                for gender in genders:
                                                        for sample in intgen_sampdict[sup_pop][population][gender]:
                                                                samples.add(sample)
                                                break
                                        
                                #Если популяция вообще не нашлась в панели,
                                #значит, она на нашу планету пока не заселилась:).
                                else:
                                        print(f'{population} - невалидное название 1000 Genomes-популяции')
                                        sys.exit()
                                        
        #Получаем и возвращаем индексы тех элементов шапки 1000
        #Genomes-VCF, которые являются названиями сэмплов,
        #отобранных в соответствии с пользовательским запросом.
        sample_indices = [header_row.index(sample) for sample in samples if sample in header_row]
        return sample_indices
