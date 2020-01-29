__version__ = 'V1.0'

import sys
sys.dont_write_bytecode = True
from backend.init_intgen_db import intgen_sampcoll_obj

def get_sample_names(populations, genders):
        '''
        Получение сэмплов, соответствующих
        заданным популяциям и полам.
        '''
        
        #Формируем часть запроса, отбирающую
        #сэмплы по принадлежности к полу.
        #Если исследователь желает работать
        #со всеми сэмплами, то расширять
        #запрос дополнительными правилами
        #фильтрации не понадобится.
        query = {'gender': {'$in': genders}}
        
        #Исследователь может по ошибке
        #указать популяцию дважды или
        #сразу суперпопуляцию и
        #принадлежащую ей субпопуляцию.
        #Конструкция с двойным ИЛИ ($in,
        #по своей сути, тоже ИЛИ) в
        #обоих случаях спасёт от дублей
        #в конечном наборе имён сэмплов.
        #Во втором случае суперпопуляция
        #как бы поглотит субпопуляцию.
        if populations != ['ALL']:
                query['$or'] = [{'super_pop': {'$in': populations}},
                                {'pop': {'$in': populations}}]
                
        #Отбор имён сэмплов.
        sample_names = [doc['sample'] for doc in intgen_sampcoll_obj.find(query)]
        
        return sample_names
