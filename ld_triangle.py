__version__ = 'V6.1'

print('''
Программа, строящая LD-матрицы для всех пар каждого
набора SNP в виде треугольной тепловой карты и/или таблицы.

Автор: Платон Быкадоров (platon.work@gmail.com), 2018-2020.
Версия: V6.1.
Лицензия: GNU General Public License version 3.
Поддержать проект: https://money.yandex.ru/to/41001832285976
Документация: https://github.com/PlatonB/ld-tools/blob/master/README.md

Обязательно!
Перед запуском программы нужно установить
MongoDB и PyMongo (см. README).

Поддерживаемые исходные файлы - таблицы,
содержащие столбец с набором refSNPIDs.
Если таких столбцов - несколько,
программа будет использовать самый левый.

Если настройки, запрашиваемые в рамках интерактивного
диалога, вам непонятны - пишите, пожалуйста, в Issues.
''')

def check_input(var):
        '''
        Проверка, правильно ли исследователь ответил на
        запрос, требующий подтверждения или отрицания.
        В случае ошибки работа программы завершится.
        '''
        if var not in ['yes', 'y', 'no', 'n', '']:
                print(f'{var} - недопустимая опция')
                sys.exit()
                
####################################################################################################

import sys, os, re, copy

#Подавление формирования питоновского кэша с
#целью предотвращения искажения результатов.
sys.dont_write_bytecode = True

from backend.init_intgen_db import *
from backend.create_intgen_db import create_intgen_db
from backend.samp_manager import get_sample_names, get_phased_genotypes
from backend.ld_calc import ld_calc
import plotly.graph_objs as go, plotly.figure_factory as ff

src_dir_path = input('\nПуть к папке с исходными файлами: ')

trg_top_dir_path = input('\nПуть к папке для конечных файлов: ')

output = input('''\nВ каком виде выводить матрицы значений
неравновесия по сцеплению (далее - LD-матрицы)?
[table(|t)|heatmap(|h)|both(|<enter>)]: ''')

if output in ['heatmap', 'h', 'both', '']:
        texts = input('''\nВыводить на диаграмму текстовую информацию?
(не рекомендуется, если строите матрицу > ~50x50 элементов; в любом
случае данные будут появляться по наведению курсора на квадратик)
[yes(|y)|no(|n|<enter>)]: ''')
        check_input(texts)
        
        if texts in ['yes', 'y']:
                val_font_size = input('''\nРазмер шрифта значений LD в квадратиках
[...|11|12(|default|<enter>)|13|...]: ''')
                if val_font_size == '':
                        val_font_size = 'default'
                elif val_font_size != 'default':
                        val_font_size = int(val_font_size)
                        
                lab_font_size = input('''\nРазмер шрифта подписей к осям
[...|11|12(|default|<enter>)|13|...]: ''')
                if lab_font_size == '':
                        lab_font_size = 'default'
                elif lab_font_size != 'default':
                        lab_font_size = int(lab_font_size)
                        
        color_pal = input('''\nЦветовая палитра тепловой карты
(https://github.com/PlatonB/ld-tools#подходящие-цветовые-палитры)
[greens(|<enter>)|amp|tempo|brwnyl|turbid|...]: ''')
        if color_pal == '':
                color_pal = 'greens'
                
elif output not in ['table', 't']:
        print(f'{output} - недопустимая опция')
        sys.exit()
        
num_of_headers = input('''\nКоличество не обрабатываемых строк
в начале каждой исходной таблицы
(игнорирование ввода ==> хэдеров/шапок в таблицах нет)
[0(|<enter>)|1|2|...]: ''')
if num_of_headers == '':
        num_of_headers = 0
else:
        num_of_headers = int(num_of_headers)
        
populations = input('''\nДля индивидов каких супер-/субпопуляций считать LD?
(http://www.internationalgenome.org/faq/which-populations-are-part-your-study/)
(несколько - через пробел)
[all(|<enter>)|eur|jpt|jpt amr yri|...]: ''').upper().split()
if populations == []:
        populations = ['ALL']
        
genders = input('''\nДля индивидов какого пола считать LD?
[male(|m)|female(|f)|both(|<enter>)]: ''').split()
if genders == ['m']:
        genders = ['male']
elif genders == ['f']:
        genders = ['female']
elif genders in [[], ['both']]:
        genders = ['male', 'female']
elif genders not in [['male'], ['female']]:
        print(f'{" ".join(genders)} - недопустимая опция')
        sys.exit()
        
ld_filter = input('''\nОбнулять значения LD, если
они меньше определённого порога?
[yes(|y)|no(|n|<enter>)]: ''')
check_input(ld_filter)

if ld_filter in ['yes', 'y']:
        thres_ld_measure = input('''\nМера LD для выставления порога
[r_square(|r)|d_prime(|d)]: ''')
        if thres_ld_measure == 'r':
                thres_ld_measure = 'r_square'
        elif thres_ld_measure == 'd':
                thres_ld_measure = 'd_prime'
        elif thres_ld_measure not in ['r_square', 'd_prime']:
                print(f'{thres_ld_measure} - недопустимая опция')
                sys.exit()
                
        ld_value_thres = float(input(f'\n{thres_ld_measure} >= '))
        
ld_measure = input('''\nМера LD для построения матриц
[r_square(|r)|d_prime(|d)]: ''')
if ld_measure == 'r':
        ld_measure = 'r_square'
elif ld_measure == 'd':
        ld_measure = 'd_prime'
elif ld_measure not in ['r_square', 'd_prime']:
        print(f'{ld_measure} - недопустимая опция')
        sys.exit()
        
#Отдельным модулем ранее были созданы
#объекты базы данных и коллекции сэмплов.
#Но существование объектов не означает, что
#база уже заполнена данными 1000 Genomes.
#Проверяем это по наличию созданных коллекций.
#Если нет ни одной коллекции, то запускаем
#модуль построения 1000 Genomes-базы с нуля.
if intgen_db_obj.command('dbstats')['collections'] == 0:
        create_intgen_db()
        
#Вызов функции, производящей отбор
#сэмплов, относящихся к указанным
#исследователем популяциям и полам.
sample_names = get_sample_names(populations,
                                genders)

#Получение списка имён всех коллекций.
intgen_coll_names = intgen_db_obj.list_collection_names()

##Работа с исходными файлами, создание конечных папок.

src_file_names = os.listdir(src_dir_path)
for src_file_name in src_file_names:
        src_file_base = '.'.join(src_file_name.split('.')[:-1])
        
        #Открытие файла исследователя на чтение.
        with open(os.path.join(src_dir_path, src_file_name)) as src_file_opened:
                
                #Считываем строки-хэдеры, чтобы сместить
                #курсор к началу основной части таблицы.
                for header_index in range(num_of_headers):
                        src_file_opened.readline()
                        
                print(f'''\n~
\n{src_file_name}
\nПодготовка набора refSNPIDs''')
                
                #refSNPIDs из обрабатываемого файла
                #будут накапливаться во множество,
                #чтобы не допустить повторяющихся.
                #С группами одинаковых снипов
                #матрицы выглядели бы абсурдно.
                rs_ids = set()
                
                #Построчное прочтение основной
                #части исходной таблицы.
                for line in src_file_opened:
                        
                        #Попытка получения refSNPID из текущей строки.
                        try:
                                rs_id = re.search(r'rs\d+', line).group()
                        except AttributeError:
                                continue
                        
                        #Добавление идентификатора во множество.
                        rs_ids.add(rs_id)
                        
        #Множества невозможно использовать
        #в запросах к MongoDB-коллекциям.
        #Поэтому множество идентификаторов
        #придётся сконвертировать в список.
        rs_ids = list(rs_ids)
        
        #Если снипов менее 2 штук, то как
        #вы себе матрицу представите?:)
        if len(rs_ids) < 2:
                print('\tМенее двух refSNPIDs')
                continue
        
        #Пользовательские refSNPID далее
        #будут разбиты по хромосомам.
        #В этот словарь будут накапливаться
        #списки полученных из базы словарей,
        #распределяясь по его разным ключам на
        #основании принадлежности SNP хромосомам.
        data_by_chrs = {}
        
        #База была спроектирована так,
        #чтобы в каждой коллекции
        #были данные одной хромосомы.
        #Коллекцию сэмплов я здесь,
        #разумеется, не имею в виду.
        #Набор refSNPID будет
        #искаться во всех хромосомных
        #коллекциях по очереди.
        #Оператор $in как бы
        #встраивает ИЛИ между
        #элементами искомого
        #списка, позволяя выдать
        #вхождения лишь его части.
        for intgen_coll_name in intgen_coll_names:
                if intgen_coll_name == 'samples':
                        continue
                elif intgen_db_obj[intgen_coll_name].count_documents({'ID': {'$in': rs_ids}}) >= 2:
                        data_by_chrs[intgen_coll_name] = list(intgen_db_obj[intgen_coll_name].find({'ID':
                                                                                                    {'$in':
                                                                                                     rs_ids}}))
                        
        #Если в коллекциях не нашлось ни одного
        #идентификатора, значит, с качеством текущего
        #набора данных исследователю совсем не повезло:(.
        #Словарь, предназначенный для похромосомного
        #распределения данных, останется пустым.
        if data_by_chrs == {}:
                print(f'\tНет валидных refSNPID')
                continue
        
        #2 или больше refSNPID накопилось, значит,
        #есть смысл создавать папку, в которую
        #будут сохраняться разделённые по хромосомам
        #текстовые и/или графические LD-матрицы.
        #Конструируем путь к этой папке и создаём её.
        trg_dir_path = os.path.join(trg_top_dir_path, f'{src_file_base}_LD_matr')
        os.mkdir(trg_dir_path)
        
        #Внутри этого цикла работа будет проводиться
        #для данных каждой хромосомы по-отдельности.
        for chr_name in data_by_chrs:
                
                #Чтобы потом проще было визуально
                #оценивать влияние физического
                #расстояния на LD, refSNPIDs надо
                #отсортировать по геномным позициям.
                #Сортируем список словарей текущей
                #хромосомы по значениям ключа POS.
                #Итоговый порядок элементов сильно
                #зависит от типа сортируемых данных.
                #Позиции предусмотрительно были
                #сконвертированы в int ещё во
                #время формирования базы,
                #поэтому на данном этапе
                #это делать уже не нужно.
                data_by_chrs[chr_name].sort(key=lambda snp_data: snp_data['POS'])
                
                #Получение отсортированных списков позиций и
                #идущих в соответствующем порядке refSNPIDs.
                poss_srtd = [snp_data['POS'] for snp_data in data_by_chrs[chr_name]]
                rs_ids_srtd = [snp_data['ID'] for snp_data in data_by_chrs[chr_name]]
                
                #Определение количества refSNPIDs, а,
                #значит, и refSNPID-содержащих словарей.
                snps_quan = len(rs_ids_srtd)
                
                print(f'\n{chr_name}: формирование LD-матрицы')
                
                #Основой текстовой или
                #графической матрицы
                #будет двумерный массив
                #такой структуры:
                '''
                 0    0    0    ...
                val   0    0    ...
                val  val   0    ...
                ...  ...  ...   ...
                '''
                
                #Построение шаблона двумерного массива
                #размером n x n, состоящего из нулей.
                #Нули в дальнейшем могут заменяться на значения LD.
                ld_two_dim = [[0 for cell_index in range(snps_quan)] for row_index in range(snps_quan)]
                
                #В случае, если будет строиться диаграмма,
                #такой же шаблон понадобится для создания
                #матрицы сопутствующей информации.
                if 'color_pal' in locals():
                        info_two_dim = copy.deepcopy(ld_two_dim)
                
                #Перебор чисел, которые будут представлять
                #собой индексы строк двумерного массива.
                for row_index in range(snps_quan):
                        
                        #Перебор чисел, которые будут служить
                        #индексами столбцов двумерного массива.
                        for col_index in range(snps_quan):
                                
                                #Матрица, в принципе, может
                                #быть квадратом, состоящим
                                #из двух одинаковых по форме
                                #и содержимому прямоугольных
                                #треугольников, разделённых
                                #диагональю 0-ячеек.
                                #Думаю, разумнее оставить лишь
                                #один из этих треугольников.
                                #Получаем только те значения,
                                #которые соответствуют ячейкам
                                #двумерного массива, индекс строки
                                #которых больше индекса столбца.
                                if row_index > col_index:
                                        
                                        #Отбор фазированных генотипов по сэмплам.
                                        snp_1_phased_genotypes = get_phased_genotypes(data_by_chrs[chr_name][row_index],
                                                                                      sample_names)
                                        snp_2_phased_genotypes = get_phased_genotypes(data_by_chrs[chr_name][col_index],
                                                                                      sample_names)
                                        
                                        #Обращение к оффлайн-калькулятору для
                                        #получения словаря с r2, D' и частотами
                                        #минорного аллеля текущей пары SNP.
                                        trg_vals = ld_calc(snp_1_phased_genotypes,
                                                           snp_2_phased_genotypes)
                                        
                                        #Для диаграммы каждое значение
                                        #LD аннотируется: параллельно с
                                        #накоплением массива LD-значений
                                        #растёт массив дополнительной
                                        #информации по каждой паре SNP.
                                        if 'info_two_dim' in locals():
                                                info_two_dim[row_index][col_index] = f'''
r2: {trg_vals["r_square"]}<br>
D': {trg_vals["d_prime"]}<br>
abs_dist: {abs(poss_srtd[col_index] - poss_srtd[row_index])}<br><br>
x.hg38_pos: {poss_srtd[col_index]}<br>
y.hg38_pos: {poss_srtd[row_index]}<br><br>
x.rsID: {rs_ids_srtd[col_index]}<br>
y.rsID: {rs_ids_srtd[row_index]}<br><br>
x.alt_freq: {trg_vals['snp_2_alt_freq']}<br>
y.alt_freq: {trg_vals['snp_1_alt_freq']}
'''
                                                
                                        #Исследователь мог установить нижний порог LD.
                                        #Соответствующий блок кода неспроста расположен после
                                        #блока накопления аннотаций: на диаграммах клеточки с
                                        #подпороговыми LD будут закрашены как нулевые, но
                                        #зато при наведении курсора там отобразятся настоящие
                                        #LD-значения, как раз извлекаемые из массива с аннотациями.
                                        #При обратном расположении этих блоков аннотации подпороговых
                                        #LD не сохранялись бы, ведь в блоке фильтрации - continue.
                                        if 'thres_ld_measure' in locals():
                                                if trg_vals[thres_ld_measure] < ld_value_thres:
                                                        continue
                                                
                                        #Если значение LD не отсеилось как
                                        #подпороговое, то попадёт в LD-матрицу:
                                        #0-ячейка будет заменена на найденное
                                        #значение LD выбранной величины.
                                        ld_two_dim[row_index][col_index] = trg_vals[ld_measure]
                                        
                #Исследователь выбрал опцию создавать
                #текстовые версии LD-матриц.
                if output in ['table', 't', 'both', '']:
                        
                        print(f'{chr_name}: сохранение текстовой LD-матрицы')
                        
                        #Создание текстового конечного файла.
                        #Прописываем в него хэдер с общими
                        #характеристиками таблицы, пустую
                        #строку и две шапки: одна - с
                        #refSNPIDs, другая - с позициями.
                        #Потом прописываем LD-строки,
                        #добавляя перед каждой из
                        #них тоже refSNPID и позицию.
                        tsv_file_name = f'{src_file_base}_{chr_name}_{ld_measure[0]}_tab.tsv'
                        with open(os.path.join(trg_dir_path, tsv_file_name), 'w') as tsv_file_opened:
                                tab, poss_srtd = '\t', list(map(lambda pos: str(pos), poss_srtd))
                                tsv_file_opened.write(f'##General\tinfo:\t{ld_measure}\t{chr_name}\t{tab.join(populations)}\t{tab.join(genders)}\n\n')
                                tsv_file_opened.write('RefSNPIDs\t\t' + '\t'.join(rs_ids_srtd) + '\n')
                                tsv_file_opened.write('\tPositions\t' + '\t'.join(poss_srtd) + '\n')
                                for row_index in range(snps_quan):
                                        line = '\t'.join([str(cell) for cell in ld_two_dim[row_index]]) + '\n'
                                        tsv_file_opened.write(rs_ids_srtd[row_index] + '\t' +
                                                              poss_srtd[row_index] + '\t' +
                                                              line)
                                        
                #Если переменная, ссылающаяся на значение, которое
                #должно стать одним из аргументов функции построения
                #диаграммы, существует, то это означает, что
                #исследователь предпочёл LD-матрицы визуализировать.
                if 'color_pal' in locals():
                        
                        print(f'{chr_name}: визуализация LD-матрицы')
                        
                        #Plotly позволяет строить тепловые карты
                        #как без надписей, так и с таковыми.
                        #Под надписями в рамках нашей задачи
                        #подразумеваются значения LD в квадратиках
                        #тепловой карты и refSNPID-лейблы осей X, Y.
                        #Допустим, исследователь предпочёл не выводить
                        #на диаграмму никаких текстовых данных.
                        #Тогда нужно будет создать обычную,
                        #не аннотированную, тепловую карту
                        #и подавить размещение лейблов осей.
                        #Для построения объекта тепловой карты,
                        #в квадратиках которой не будет текста,
                        #программа использует класс Heatmap.
                        #Класс принимает посчитанные LD и другую
                        #информацию по парам SNPs, заданную
                        #исследователем цветовую схему и два
                        #константных аргумента: наличие расстояния
                        #между квадратиками, отсутствие цветовой шкалы.
                        #Объект диаграммы дополняется словареподобным
                        #объектом класса Layout с настройками
                        #осей: для начала - запретом вывода лейблов.
                        if texts in ['no', 'n', '']:
                                trace = go.Heatmap(z=ld_two_dim,
                                                   hovertext=info_two_dim,
                                                   hoverinfo='text',
                                                   xgap=1,
                                                   ygap=1,
                                                   colorscale=color_pal,
                                                   showscale=False)
                                layout = go.Layout(xaxis={'showticklabels': False},
                                                   yaxis={'showticklabels': False})
                                ld_heatmap = go.Figure(data=trace,
                                                       layout=layout)
                                
                        #Исследователь дал добро выводить на диаграмму надписи.
                        else:
                                
                                #Для создания тепловых карт со значениями в
                                #квадратиках (аннотированных) Plotly предоставляет
                                #высокоуровневую функцию create_annotated_heatmap.
                                #Объект аннотированной тепловой карты представляет
                                #собой структуру со свойствами словаря, в которую
                                #вложены структуры со свойствами словарей и
                                #списков, содержащие параметры диаграммы.
                                #С помощью аргументов create_annotated_heatmap
                                #можно менять настройки диаграммы, находящиеся
                                #в значении ключа верхнего уровня 'data'.
                                #Самые основные аргументы этой функции - массивы
                                #с LD-данными, refSNPIDs в качестве лейблов по X, Y
                                #и дополнительной информацией по каждой паре SNPs.
                                #Аргумент, который определяется на этапе интерактивного
                                #диалога исследователем - цветовая схема тепловой карты.
                                #Не изменяемый в рамках диалога аргумент - наличие
                                #разделительных линий между квадратиками.
                                #В create_annotated_heatmap по умолчанию
                                #не допускается вывод цветовой шкалы.
                                #Так и оставим ради экономии пространства.
                                ld_heatmap = ff.create_annotated_heatmap(ld_two_dim,
                                                                         x=rs_ids_srtd,
                                                                         y=rs_ids_srtd,
                                                                         hovertext=info_two_dim,
                                                                         hoverinfo='text',
                                                                         xgap=1,
                                                                         ygap=1,
                                                                         colorscale=color_pal)
                                
                                #Исследователь кастомизировал размер
                                #шрифта значений внутри квадратиков.
                                #Тогда в каждый подсловарь структуры
                                #'annotations' будет добавлена пара
                                #ключ-значение с данным размером.
                                if val_font_size != 'default':
                                        for ann_num in range(len(ld_heatmap['layout']['annotations'])):
                                                ld_heatmap['layout']['annotations'][ann_num]['font']['size'] = val_font_size
                                                
                                #При желании, исследователь мог задать недефолтное
                                #значение параметра размера шрифта лейблов осей.
                                if lab_font_size != 'default':
                                        ld_heatmap['layout']['xaxis']['tickfont'] = {'size': lab_font_size}
                                        ld_heatmap['layout']['yaxis']['tickfont'] = {'size': lab_font_size}
                                        
                        #Заголовок будет добавляться в любую
                        #диаграмму - обычную и аннотированную.
                        #Он должен содержать в себе единые для
                        #всей тепловой карты характеристики:
                        #величину LD (r2 или D'), задающуя
                        #цвета квадратиков, номер хромосомы, в
                        #которой локализуются SNP, формирующие
                        #оси диаграммы, а также названия
                        #популяций и полов, по которым ранее
                        #отбирались сэмплы для рассчёта LD.
                        ld_heatmap['layout']['title'] = {'text':
f'''defines_color: {ld_measure} ░
chrom: {chr_name} ░
populations: {", ".join(populations)} ░
genders: {", ".join(genders)}'''}
                        
                        #В тепловых картах Plotly, обычных и аннотированных,
                        #ось Y по умолчанию направлена снизу вверх.
                        #Чтобы ориентация тепловых карт, создаваемых
                        #этой программой, была такая же, как в
                        #LDMatrix, переворачиваем диаграмму по Y.
                        #Параметры осей заложены в структуры,
                        #принадлежащие ключу верхнего уровня 'layout'.
                        ld_heatmap['layout']['yaxis']['autorange'] = 'reversed'
                        
                        #Построение диаграммы и её сохранение в HTML.
                        html_file_path = f'{os.path.join(trg_dir_path, src_file_base)}_{chr_name}_{ld_measure[0]}_diag.html'
                        ld_heatmap.write_html(html_file_path)
                        
client.close()
