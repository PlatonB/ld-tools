__version__ = 'V1.0'

import sqlite3

def get_sample_names(gend_names, pop_names, intgen_convdb_path):
        '''
        Получение сэмплов, соответствующих
        заданным популяциям и полам.
        '''
        
        #Формируем часть запроса, отбирающую
        #сэмплы по принадлежности к полу.
        #Если исследователь желает работать
        #со всеми сэмплами, то далее
        #расширять запрос дополнительными
        #правилами фильтрации не понадобится.
        query = f'SELECT sample FROM samples WHERE gender IN {gend_names}'
        
        #Исследователь может по ошибке
        #указать популяцию дважды или
        #сразу суперпопуляцию и
        #принадлежащую ей субпопуляцию.
        #Конструкция с двойным ИЛИ (IN,
        #по своей сути, тоже ИЛИ) в
        #обоих случаях спасёт от дублей
        #в конечном наборе имён сэмплов.
        #Во втором случае суперпопуляция
        #как бы поглотит субпопуляцию.
        if pop_names != ('ALL',):
                query += f' AND (super_pop IN {pop_names} OR pop IN {pop_names})'
        query = query.replace(',)', ')')
        
        #Инициализация базы данных, одна
        #из таблиц которой содержит сэмплы,
        #распределённые по полам и популяциям.
        with sqlite3.connect(intgen_convdb_path) as conn:
                cursor = conn.cursor()
                
                #Отбор имён сэмплов.
                sample_names = [sample_tup[0] for sample_tup in cursor.execute(query)]
                
                #Дисконнектимся.
                cursor.close()
                
        return sample_names
