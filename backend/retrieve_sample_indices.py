__version__ = 'V1.1'

import json, dbm, gzip, sys

def retrieve_sample_indices(intgen_sampjson_path, populations, genders, any_intgen_vcfdb_path):
        '''
        Получение сэмплов, соответствующих заданным популяциям и полам.
        '''
        
        #Открытие JSON-варианта панели сэмплов на чтение.
        #Извлечение хранящегося там объекта в оперативную память.
        with open(intgen_sampjson_path) as intgen_sampjson_opened:
                intgen_sampdict = json.load(intgen_sampjson_opened)
                
        #Во всех VCF-таблицах 1000
        #Genomes одна и та же шапка.
        #Открываем любой VCF-содержащий архив
        #на чтение и достаём из таблицы шапку.
        #Разархивируем её, конвертируем
        #в Юникод и преобразуем в список.
        with dbm.open(any_intgen_vcfdb_path) as any_intgen_vcfdb_opened:
                header_row = gzip.decompress(any_intgen_vcfdb_opened['header_line']).decode('utf-8').split('\t')
                
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
