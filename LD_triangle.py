__version__ = 'V5.2'

print('''
Программа, строящая LD-матрицы для всех пар каждого
набора SNP в виде треугольной тепловой карты и/или таблицы.

Автор: Платон Быкадоров (platon.work@gmail.com), 2018-2019.
Версия: V5.2.
Лицензия: GNU General Public License version 3.
Поддержать проект: https://money.yandex.ru/to/41001832285976
Документация: https://github.com/PlatonB/ld-tools/blob/master/README.md

Обязательно!
Перед запуском программы нужно установить модули:
sudo pip3 install plyvel plotly numpy

Поддерживаемые исходные файлы - таблицы,
содержащие столбец с набором refSNPIDs.
Если таких столбцов - несколько,
программа будет использовать самый левый.

Если настройки, запрашиваемые в рамках интерактивного
диалога, вам непонятны - пишите, пожалуйста, в Issues.
''')

def check_input(var):
        '''
        Проверка, правильно ли пользователь ответил на
        запрос, требующий подтверждения или отрицания.
        В случае ошибки работа программы завершится.
        '''
        if var not in ['yes', 'y', 'no', 'n', '']:
                print(f'{var} - недопустимая опция')
                sys.exit()
                
####################################################################################################

print('\nИмпорт модулей программы...')

import sys

#Раньше при нескольких подряд
#запусках этой программы можно было
#наблюдать, по меньшей мере, 2 бага.
#1. Чаще: новая диаграмма рисовалась
#тех же цветов, что и предыдущая.
#2. Реже, но с разрушительными последствиями:
#выдавались неправильные значения LD.
#Всё это лечится подавлением
#формирования питоновского кэша.
sys.dont_write_bytecode = True

import os, re, gzip, plyvel, copy
from backend.prepare_intgen_data import process_intgen_data
from backend.retrieve_sample_indices import retrieve_sample_indices
from backend.ld_calc import ld_calc
import plotly as py, plotly.graph_objs as go, plotly.figure_factory as ff

src_dir_path = input('\nПуть к папке с исходными файлами: ')

trg_top_dir_path = input('\nПуть к папке для конечных файлов: ')

intgen_dir_path = input('\nПуть к папке для данных 1000 Genomes: ')

output = input('''\nВ каком виде выводить матрицы значений
неравновесия по сцеплению (далее - LD-матрицы)?
(игнорирование ввода ==> в обоих)
[table(|t)|heatmap(|h)|both(|<enter>)]: ''')

if output in ['heatmap', 'h', 'both', '']:
        texts = input('''\nВыводить на диаграмму текстовую информацию?
(не рекомендуется, если строите матрицу > ~50x50 элементов; в любом
случае данные будут появляться по наведению курсора на квадратик)
(игнорирование ввода ==> не выводить)
[yes(|y)|no(|n|<enter>)]: ''')
        check_input(texts)
        
        if texts in ['yes', 'y']:
                val_font_size = input('''\nРазмер шрифта значений LD в квадратиках
(игнорирование ввода ==> значение по умолчанию)
[default(|<enter>)|...|11|12|13|...]: ''')
                if val_font_size == '':
                        val_font_size = 'default'
                elif val_font_size != 'default':
                        val_font_size = int(val_font_size)
                        
                lab_font_size = input('''\nРазмер шрифта подписей к осям
(игнорирование ввода ==> значение по умолчанию)
[default(|<enter>)|...|11|12|13|...]: ''')
                if lab_font_size == '':
                        lab_font_size = 'default'
                elif lab_font_size != 'default':
                        lab_font_size = int(lab_font_size)
                        
        color_map = input('''\nЦветовая палитра тепловой карты
(https://github.com/plotly/plotly.py/blob/bdafae81e747dfc11ea7c6c64fa9d09117aa1816/plotly/colors.py#L90)
(игнорирование ввода ==> зелёные оттенки)
[greens(|<enter>)|hot|earth|electric|...]: ''').capitalize()
        if color_map == '':
                color_map = 'Greens'
                
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
        
populations = input('''\nДля индивидов какой(-их) супер-/субпопуляции(-ий) считать LD?
(http://www.internationalgenome.org/faq/which-populations-are-part-your-study/)
(несколько - через пробел)
(игнорирование ввода ==> всех)
[all(|<enter>)|eur|jpt|jpt amr yri|...]: ''').upper().split()
if populations == []:
        populations = ['ALL']
        
genders = input('''\nДля индивидов какого пола считать LD?
(игнорирование ввода ==> обоих)
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
(игнорирование ввода ==> не
фильтровать пары по порогу LD)
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
        
#Вызов функции, которая, во-первых, скачает
#заархивированные VCF и панель сэмплов проекта
#1000 Genomes, если они ещё не скачаны, во-вторых,
#сконвертирует панель в JSON, а данные изо
#всех VCF разместит в единую LevelDB-базу.
#LevelDB-база данных (далее - база) будет
#состоять из ключей - ID снипов всех хромосом,
#и значений - сжатых строк упомянутых VCF.
intgen_sampjson_path, intgen_vcfdb_path = process_intgen_data(intgen_dir_path)[:2]

#Вызов функции, производящей отбор
#сэмплов, относящихся к указанным
#пользователем популяциям и полам,
#и возвращающей индексы этих сэмплов.
sample_indices = retrieve_sample_indices(intgen_sampjson_path,
                                         populations,
                                         genders,
                                         intgen_vcfdb_path)

#Открытие базы.
intgen_vcfdb_opened = plyvel.DB(intgen_vcfdb_path, create_if_missing=False)

##Работа с исходными файлами, создание конечных папок.

src_file_names = os.listdir(src_dir_path)
for src_file_name in src_file_names:
        src_file_base = '.'.join(src_file_name.split('.')[:-1])
        
        #Открытие файла пользователя на чтение.
        with open(os.path.join(src_dir_path, src_file_name)) as src_file_opened:
                
                #Считываем строки-хэдеры, чтобы сместить
                #курсор к началу основной части таблицы.
                for header_index in range(num_of_headers):
                        src_file_opened.readline()
                        
                print(f'''\n~
\n{src_file_name}
\nПодготовка набора SNPs. Предупреждения (если будут):''')
                
                ##Очистка пользовательского набора refSNPID от повторяющихся.
                
                #Для этого данные refSNPID будут накапливаться во множество.
                rs_ids = set()
                
                #Построчное прочтение основной части исходной таблицы.
                for line in src_file_opened:
                        
                        #Попытка получения refSNPID из текущей строки.
                        try:
                                rs_id = re.search(r'rs\d+', line).group()
                        except AttributeError:
                                continue
                        
                        #Добавление идентификатора во множество.
                        rs_ids.add(rs_id)
                        
        #Проверка, есть ли в пользовательском наборе хотя бы 2 refSNPID.
        if len(rs_ids) < 2:
                print('\nФайл содержит менее двух refSNPIDs. Будет проигнорирован')
                continue
        
        #2 или больше refSNPID накопились, значит,
        #есть смысл создавать папку, в которую
        #будут сохраняться разделённые по хромосомам
        #текстовые и/или графические LD-матрицы.
        #Конструируем путь к этой папке и создаём её.
        trg_dir_path = os.path.join(trg_top_dir_path, f'{src_file_base}_LD_matr')
        os.mkdir(trg_dir_path)
        
        ##Распределение по хромосомам пользовательских refSNPIDs.
        
        #В этот словарь будут накапливаться преобразованные в списки
        #refSNPID-содержащие строки из базы, распределяясь по разным
        #элементам словаря на основе принадлежности SNP хромосомам.
        rows_by_chrs = {}
        
        #Перебор избавленных от повторов пользовательских refSNPID.
        for rs_id in rs_ids:
                
                #Попытка извлечения из базы сжатой строки,
                #соответствующей текущему пользовательскому refSNPID.
                #Декомпрессия результата, превращение из
                #байт-строки в обычную и преобразование в список.
                #Если refSNPID в базе нет, то вместо сжатой строки вернётся None.
                #В результате разархивации None выйдет пустая байтовая строка.
                #После её конвертаций, упомянутых выше, сформируется список
                #с обычной пустой строкой в качестве единственного элемента.
                row = gzip.decompress(intgen_vcfdb_opened.get(rs_id.encode())).decode('utf-8').split('\t')
                
                #Если refSNPID из файла
                #пользователя не обнаружился
                #в базе, значит он - невалидный.
                #В матрицу он допущен не будет.
                #Одна из вероятных причин - в том,
                #что именуемый им SNP - не биаллельный.
                #Переходим к следующему refSNPID.
                if row == ['']:
                        print(f'\t{rs_id} - невалидный refSNPID (возможно, ID мультиаллельного SNP). В матрицу не пойдёт')
                        continue

                #refSNPID в базе нашёлся.
                #Извлечение из refSNPID-содержащего списка
                #номера хромосомы соответствующего SNP.
                chr_num = row[0]
                
                #Ключами ранее созданного словаря
                #будут номера хромосом, а значениями -
                #двумерные массивы, каждый элемент
                #которых - refSNPID-содержащий список.
                #Если очередной такой список относится
                #к хромосоме, которая ранее не
                #добавлялась, то создаётся новая пара
                #ключ-значение (хромосома-пустой список).
                if chr_num not in rows_by_chrs:
                        rows_by_chrs[chr_num] = []
                        
                #Размещение текущего refSNPID-содержащего списка
                #в словарь с привязкой к хромосомному ключу.
                rows_by_chrs[chr_num].append(row)
                
        #Признаком того, что невалидны все refSNPID исходного
        #файла, будет пустота словаря, предназначенного
        #для похромосомного распределения данных.
        if rows_by_chrs == {}:
                print(f'\t{src_file_name} не содержит ни одного валидного refSNPID')
                continue
                
        ##Сортировка списка refSNPID по их позиции в хромосоме.
        
        #Внутри этого цикла работа будет проводиться
        #для данных каждой хромосомы по-отдельности.
        for chr_num in rows_by_chrs:
                
                #Определение количества refSNPID-содержащих списков.
                rows_in_chr_amount = len(rows_by_chrs[chr_num])
                
                #Поскольку LD определяется для пар SNP, то если на
                #данную хромосому пришёлся только 1 SNP пользовательского
                #набора, необходимо перейти к следующей хромосоме (если имеется).
                if rows_in_chr_amount < 2:
                        continue
                
                #Сортировка двумерного массива текущей хромосомы по позициям SNPs.
                rows_by_chrs[chr_num].sort(key=lambda row: int(row[1]))
                
                #Получение списков отсортированных позиций и
                #идущих в соответствующем порядке refSNPIDs.
                poss_srtd = [row[1] for row in rows_by_chrs[chr_num]]
                rs_ids_srtd = [row[2] for row in rows_by_chrs[chr_num]]
                
                print(f'\nхромосома {chr_num}: формирование LD-матрицы...')
                
                ##Создание двумерных массивов такой структуры:
                '''
                 0      0     0
                val     0     0
                val    val    0
                '''
                
                #Построение шаблона двумерного массива
                #размером n x n, состоящего из нулей.
                #Нули в дальнейшем могут заменяться на значения LD.
                ld_two_dim = [[0 for cell_index in range(rows_in_chr_amount)] for row_index in range(rows_in_chr_amount)]
                
                #В случае, если будет строиться диаграмма,
                #такой же шаблон понадобится для создания
                #матрицы сопутствующей информации.
                if 'color_map' in locals():
                        info_two_dim = copy.deepcopy(ld_two_dim)
                
                #Перебор чисел, которые будут представлять
                #собой индексы строк двумерных массивов.
                for row_index in range(rows_in_chr_amount):
                        
                        #Перебор чисел, которые будут служить
                        #индексами столбцов двумерных массивов.
                        for col_index in range(rows_in_chr_amount):
                                
                                #Матрица, в принципе, может быть квадратом,
                                #состоящим из двух одинаковых по форме и содержимому
                                #прямоугольных треугольников, разделённых диагональю 0-ячеек.
                                #Думаю, разумнее оставить лишь один из этих треугольников.
                                #Получаем только те значения, которые соответствуют
                                #ячейкам, индекс строки которых больше индекса столбца.
                                if row_index > col_index:
                                        
                                        #Обращение к оффлайн-калькулятору для
                                        #получения словаря с r2, D' и частотами
                                        #минорного аллеля текущей пары SNP.
                                        trg_vals = ld_calc(rows_by_chrs[chr_num][row_index],
                                                           rows_by_chrs[chr_num][col_index],
                                                           sample_indices)
                                        
                                        #Для диаграммы каждое значение LD аннотируется:
                                        #параллельно с накоплением массива LD-значений растёт
                                        #массив дополнительной информации по каждой паре SNP.
                                        if 'info_two_dim' in locals():
                                                info_two_dim[row_index][col_index] = f'''
r2: {trg_vals["r_square"]}<br>
D': {trg_vals["d_prime"]}<br>
abs_dist: {abs(int(poss_srtd[col_index]) - int(poss_srtd[row_index]))}<br><br>
chr: {chr_num}<br>
x.hg38_pos: {poss_srtd[col_index]}<br>
y.hg38_pos: {poss_srtd[row_index]}<br><br>
x.rsID: {rs_ids_srtd[col_index]}<br>
y.rsID: {rs_ids_srtd[row_index]}<br><br>
x.alt_freq: {trg_vals['snp_2_alt_freq']}<br>
y.alt_freq: {trg_vals['snp_1_alt_freq']}<br><br>
pops: {", ".join(populations)}<br>
gends: {", ".join(genders)}
'''
                                                
                                        #Пользователь мог установить нижний порог LD.
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
                                        
                ##Сохранение результатов.
                
                #Пользователь указал создавать текстовые версии LD-матриц.
                if output in ['table', 't', 'both', '']:
                        
                        print(f'хромосома {chr_num}: сохранение текстовой LD-матрицы...')
                        
                        #Создание текстового конечного файла.
                        #Прописываем в него хэдер с общими
                        #характеристиками таблицы, пустую
                        #строку и две шапки: одна - с
                        #refSNPIDs, другая - с позициями.
                        #Потом прописываем LD-строки,
                        #добавляя перед каждой из
                        #них тоже refSNPID и позицию.
                        tsv_file_name = f'{src_file_base}_chr{chr_num}_{ld_measure[0]}_tab.tsv'
                        with open(os.path.join(trg_dir_path, tsv_file_name), 'w') as tsv_file_opened:
                                tab = '\t'
                                tsv_file_opened.write(f'#General\tinfo:\t{ld_measure}\tchr{chr_num}\t{tab.join(populations)}\t{tab.join(genders)}\n\n')
                                tsv_file_opened.write('RefSNPIDs\t\t' + '\t'.join(rs_ids_srtd) + '\n')
                                tsv_file_opened.write('\tPositions\t' + '\t'.join(poss_srtd) + '\n')
                                rs_id_index = 0
                                for row in ld_two_dim:
                                        line = '\t'.join([str(cell) for cell in row]) + '\n'
                                        tsv_file_opened.write(rs_ids_srtd[rs_id_index] + '\t' +
                                                              poss_srtd[rs_id_index] + '\t' +
                                                              line)
                                        rs_id_index += 1
                                        
                #Если переменная, ссылающаяся на значение, которое
                #должно стать одним из аргументов функции построения
                #диаграммы, существует, то это означает, что
                #пользователь предпочёл LD-матрицы визуализировать.
                if 'color_map' in locals():
                        
                        print(f'хромосома {chr_num}: визуализация LD-матрицы...')
                        
                        #Plotly позволяет строить тепловые карты
                        #как без надписей, так и с таковыми.
                        #Под надписями в рамках нашей задачи
                        #подразумеваются значения LD в квадратиках
                        #тепловой карты и refSNPID-лейблы осей X, Y.
                        #Допустим, пользователь предпочёл не выводить
                        #на диаграмму никаких текстовых данных.
                        #Тогда нужно будет создать обычную,
                        #не аннотированную, тепловую карту
                        #и подавить размещение лейблов осей.
                        #Для построения объекта тепловой карты,
                        #в квадратиках которой не будет текста,
                        #программа использует класс Heatmap.
                        #Класс принимает посчитанные LD и другую
                        #информацию по парам SNPs, заданную
                        #пользователем цветовую схему и ряд
                        #константных аргументов: наличие
                        #расстояния между квадратиками,
                        #потемнение клеточек по мере увеличения
                        #значений и отсутствие цветовой шкалы.
                        #Объект диаграммы дополняется словареподобным
                        #объектом класса Layout с настройками
                        #осей: для начала - запретом вывода лейблов.
                        if texts in ['no', 'n', '']:
                                trace = go.Heatmap(z=ld_two_dim,
                                                   hovertext=info_two_dim,
                                                   hoverinfo='text',
                                                   xgap=1,
                                                   ygap=1,
                                                   colorscale=color_map,
                                                   reversescale=True,
                                                   showscale=False)
                                layout = go.Layout(xaxis={'showticklabels': False},
                                                   yaxis={'showticklabels': False})
                                ld_heatmap = go.Figure(data=[trace],
                                                       layout=layout)
                                
                        #Пользователь дал добро выводить на диаграмму надписи.
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
                                #диалога пользователем - цветовая схема тепловой карты.
                                #Не изменяемые в рамках диалога аргументы: наличие
                                #разделительных линий между квадратиками, более тёмная
                                #заливка квадратиков при больших значениях в них.
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
                                                                         colorscale=color_map,
                                                                         reversescale=True)
                                
                                #Пользователь кастомизировал размер шрифта значений внутри квадратиков.
                                #Тогда в каждый подсловарь структуры 'annotations' будет
                                #добавлена пара ключ-значение с данным размером.
                                if val_font_size != 'default':
                                        for ann_num in range(len(ld_heatmap['layout']['annotations'])):
                                                ld_heatmap['layout']['annotations'][ann_num]['font']['size'] = val_font_size
                                                
                                #При желании, пользователь мог задать недефолтное
                                #значение параметра размера шрифта лейблов осей.
                                if lab_font_size != 'default':
                                        ld_heatmap['layout']['xaxis']['tickfont'] = {'size': lab_font_size}
                                        ld_heatmap['layout']['yaxis']['tickfont'] = {'size': lab_font_size}
                                        
                        #В тепловых картах Plotly, обычных и аннотированных,
                        #ось Y по умолчанию направлена снизу вверх.
                        #Чтобы ориентация тепловых карт, создаваемых
                        #этой программой, была такая же, как в
                        #LDMatrix, переворачиваем диаграмму по Y.
                        #Параметры осей заложены в структуры,
                        #принадлежащие ключу верхнего уровня 'layout'.
                        ld_heatmap['layout']['yaxis']['autorange'] = 'reversed'
                        
                        #Построение диаграммы и её сохранение в HTML.
                        html_file_path = f'{os.path.join(trg_dir_path, src_file_base)}_chr{chr_num}_{ld_measure[0]}_diag.html'
                        ld_heatmap.write_html(html_file_path)
                        
intgen_vcfdb_opened.close()
