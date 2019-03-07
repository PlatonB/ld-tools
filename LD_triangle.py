__version__ = 'V3.1'

print('''
Программа, строящая LD-матрицы для всех пар каждого
набора SNP в виде треугольной тепловой карты и/или таблицы.

Автор: Платон Быкадоров (platon.work@gmail.com), 2018-2019.
Версия: V3.1.
Лицензия: GNU General Public License version 3.
Поддержать проект: https://money.yandex.ru/to/41001832285976
Документация: https://github.com/PlatonB/ld-tools/blob/master/README.md

Обязательно! Установка библиотеки визуализации:
sudo pip install plotly

Поддерживаемые исходные файлы - таблицы, содержащие столбец с набором refSNPIDs.
Если таких столбцов - несколько, программа будет использовать самый левый.

Если настройки, запрашиваемые в рамках интерактивного
диалога, вам непонятны - пишите, пожалуйста, в Issues.
''')

def check_input(var):
        '''
        Проверка, правильно ли пользователь ответил на
        запрос, требующий подтверждения или отрицания.
        В случае ошибки работа программы завершится.
        '''
        if var != 'yes' and var != 'y' and var != 'no' \
           and var != 'n' and var != '':
                print(f'{var} - недопустимая опция')
                sys.exit()
                
####################################################################################################

print('\nИмпорт модулей программы...')

import sys, random, os, re, gzip, dbm, shutil
sys.path.insert(0, os.path.join(os.getcwd(), 'backend'))
from prepare_intgen_data import process_intgen_data
from retrieve_sample_indices import retrieve_sample_indices
from ld_calc import ld_calc
import plotly as py, plotly.figure_factory as ff

src_dir_path = input('\nПуть к папке с исходными файлами: ')
trg_top_dir_path = input('\nПуть к папке для конечных файлов: ')
intgen_dir_path = input('\nПуть к папке для данных 1000 Genomes: ')

output = input('''\nВ каком виде выводить матрицы значений
неравновесия по сцеплению (далее - LD-матрицы)?
(игнорирование ввода ==> в обоих)
[table|heatmap|both(|<enter>)]: ''')
if output != 'table' and output != 'heatmap' and output != 'both' and output != '':
        print(f'{output} - недопустимая опция')
        sys.exit()
        
elif output == 'heatmap' or output == 'both' or output == '':
        texts = input('''\nВыводить на диаграмму текстовую информацию?
(не рекомендуется, если строите матрицу > ~50x50 элементов; в любом
случае данные будут появляться по наведению курсора на квадратик)
(игнорирование ввода ==> не выводить)
[yes(|y)|no(|n|<enter>)]: ''')
        check_input(texts)
        
        if texts == 'yes' or texts == 'y':
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
[Greens(|<enter>)|Hot|Earth|Electric|...]: ''')
        if color_map == '':
                color_map = 'Greens'
                
num_of_headers = input('''\nКоличество не обрабатываемых строк
в начале каждой исходной таблицы
(игнорирование ввода ==> хэдеров/шапок в таблицах нет)
[0(|<enter>)|1|2|...]: ''')
if num_of_headers == '':
        num_of_headers = 0
else:
        num_of_headers = int(num_of_headers)
        
populations = input('''\nВ пределах какой(-их) супер-/субпопуляции(-ий) считать LD?
(http://www.internationalgenome.org/faq/which-populations-are-part-your-study/)
(несколько - через пробел)
(игнорирование ввода ==> всех)
[ALL(|<enter>)|EUR|JPT|AMR YRI|...]: ''').split()
if populations == []:
        populations = ['ALL']
        
genders = input('''\nСреди индивидов какого пола считать LD?
(игнорирование ввода ==> обоих)
[male|female|both(|<enter>)]: ''').split()
if genders == [] or genders == ['both']:
        genders = ['male', 'female']
elif genders != ['male'] and genders != ['female']:
        print(f'{genders[0]} - недопустимая опция')
        sys.exit()
        
ld_filter = input('''\nОбнулять значения LD, если
они меньше определённого порога?
(игнорирование ввода ==> не
фильтровать пары по порогу LD)
[yes(|y)|no(|n|<enter>)]: ''')
check_input(ld_filter)

if ld_filter == 'yes' or ld_filter == 'y':
        thres_ld_measure = input('''\nМера LD для выставления порога
[r_square(|r)|d_prime(|d)]: ''')
        if thres_ld_measure == 'r':
                thres_ld_measure = 'r_square'
        elif thres_ld_measure == 'd':
                thres_ld_measure = 'd_prime'
        elif thres_ld_measure != 'r_square' and thres_ld_measure != 'd_prime':
                print(f'{thres_ld_measure} - недопустимая опция')
                sys.exit()
                
        ld_value_thres = float(input(f'\n{thres_ld_measure} >= '))
        
ld_measure = input('''\nМера LD для построения матриц
[r_square(|r)|d_prime(|d)]: ''')
if ld_measure == 'r':
        ld_measure = 'r_square'
elif ld_measure == 'd':
        ld_measure = 'd_prime'
elif ld_measure != 'r_square' and ld_measure != 'd_prime':
        print(f'{ld_measure} - недопустимая опция')
        sys.exit()
        
#Вызов функции, которая скачает заархивированные VCF и панель сэмплов проекта
#1000 Genomes, если они ещё не скачаны, а также разместит строки изо всех VCF
#в единую dbm-базу и вернёт абсолютные пути ко всем скачанным и созданным файлам.
#Эта dbm-база данных (далее - база) будет состоять из ключей -
#refSNPID всех хромосом, и значений - сжатых строк упомянутых VCF.
intgen_sampjson_path, intgen_natgz_paths, intgen_natdb_paths = process_intgen_data(intgen_dir_path)

#Удаление путей к 1000 Genomes-архивам, т.к. у нас уже
#есть база, из которой легко извлечь refSNP-содержащие строки,
#а исходные архивы именно в этой программе более не пригодятся.
del intgen_natgz_paths

#Вызов функции, производящей отбор сэмплов, относящихся к указанным
#пользователем популяциям и полам, и возвращающей индексы этих сэмплов.
sample_indices = retrieve_sample_indices(intgen_sampjson_path, populations, genders, random.choice(intgen_natdb_paths))

##Работа с исходными файлами, создание конечных папок.

src_file_names = os.listdir(src_dir_path)
for src_file_name in src_file_names:
        src_file_base = '.'.join(src_file_name.split('.')[:-1])
        
        #Создание папки, в которую будут сохраняться разделённые
        #по хромосомам текстовые и/или графические LD-матрицы.
        trg_dir_path = os.path.join(trg_top_dir_path, f'{src_file_base}_LD_matr')
        os.mkdir(trg_dir_path)
        
        #Открытие файла пользователя на чтение.
        with open(os.path.join(src_dir_path, src_file_name)) as src_file_opened:
                
                #Считываем строки-хэдеры, чтобы сместить
                #курсор к началу основной части таблицы.
                for header_index in range(num_of_headers):
                        src_file_opened.readline()
                        
                print(f'''\n\n{src_file_name}: подготовка набора SNPs...
Предупреждения (если будут):''')
                
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
                        
                #Проверка, есть ли пользовательском наборе хотя бы 2 refSNPID.
                if len(rs_ids) < 2:
                        print(f'{src_file_name} содержит менее двух refSNPIDs. Будет проигнорирован')
                        continue
                        
                ##Распределение по хромосомам пользовательских refSNPIDs.
                
                #В этот словарь будут накапливаться преобразованные в списки
                #refSNPID-содержащие строки базы, распределяясь по разным
                #элементам словаря на основе принадлежности SNP хромосомам.
                rows_by_chrs = {}
                
                #Перебор избавленных от повторов пользовательских refSNPID.
                for rs_id in rs_ids:
                        
                        #Перебор баз с целью найти refSNPID-содержащую строку.
                        for intgen_natdb_path in intgen_natdb_paths:
                                
                                #Открытие каждой базы на чтение.
                                with dbm.open(intgen_natdb_path) as intgen_natdb_opened:
                                        
                                        #Попытка извлечения из очередной базы сжатой строки,
                                        #соответствующей текущему пользовательскому refSNPID.
                                        #Декомпрессия результата, превращение из
                                        #байт-строки в обычную и преобразование в список.
                                        #Если refSNPID в базе нет, то вместо сжатой строки вернётся None.
                                        #В результате разархивации None выйдет пустая байтовая строка.
                                        #После её конвертаций, упомянутых выше, сформируется список
                                        #с обычной пустой строкой в качестве единственного элемента.
                                        row = gzip.decompress(intgen_natdb_opened.get(rs_id)).decode('utf-8').split('\t')
                                        
                                        #Если содержимое списка получилось
                                        #иное, нежели пустая строка, значит
                                        #refSNPID в одной из баз нашёлся.
                                        #Прерываем поиск.
                                        if row != ['']:
                                                break
                                        
                        #Если refSNPID из файла
                        #пользователя не обнаружился
                        #в базах, значит он - невалидный.
                        #В матрицу он допущен не будет.
                        #Одна из вероятных причин - в том,
                        #что именуемый им SNP - не биаллельный.
                        else:
                                print(f'{rs_id} - невалидный refSNPID (возможно, ID мультиаллельного SNP). В матрицу не пойдёт')
                                continue
                        
                        #Извлечение из refSNPID-содержащего списка
                        #номера хромосомы соответствующего SNP.
                        chr_num = row[0]
                        
                        #Ключами ранее созданного словаря будут номера хромосом, а значениями -
                        #двумерные массивы, каждый элемент которых - refSNPID-содержащий список.
                        #Если очередной такой список относится к хромосоме, которая ранее не
                        #добавлялась, то создаётся новая пара ключ-значение (хромосома-пустой список).
                        if chr_num not in rows_by_chrs:
                                rows_by_chrs[chr_num] = []
                                
                        #Размещение текущего refSNPID-содержащего списка
                        #в словарь с привязкой к хромосомному ключу.
                        rows_by_chrs[chr_num].append(row)
                        
                #Признаком того, что невалидны все refSNPID исходного
                #файла, будет пустота словаря, предназначенного
                #для похромосомного распределения данных.
                if rows_by_chrs == {}:
                        print(f'{src_file_name} не содержит ни одного валидного refSNPID')
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
                        
                        #Получение отсортированного по позиции списка refSNPID.
                        rs_ids_srtd = [row[2] for row in rows_by_chrs[chr_num]]
                        
                        print(f'\n{src_file_name}: хромосома {chr_num}: формирование LD-матрицы...')
                        
                        ##Создание двумерного массива такой структуры:
                        '''
                          0        0      0
                        LDval      0      0
                        LDval    LDval    0
                        '''
                        
                        #Построение шаблона двумерного массива
                        #размером n x n, состоящего из нулей.
                        #Нули в дальнейшем будут заменяться на значения LD.
                        ld_two_dim = [[0 for cell_index in range(rows_in_chr_amount)] for row_index in range(rows_in_chr_amount)]
                        
                        #Перебор чисел, которые будут представлять
                        #собой индексы строк двумерного массива.
                        for row_index in range(rows_in_chr_amount):
                                
                                #Перебор чисел, которые будут служить
                                #индексами столбцов двумерного массива.
                                for col_index in range(rows_in_chr_amount):
                                        
                                        #Матрица значений может представлять собой квадрат, состоящий
                                        #из двух одинаковых по форме и содержимому прямоугольных
                                        #треугольника, разделённых диагональю 0-ячеек.
                                        #Думаю, разумнее оставить лишь один из этих треугольников.
                                        #Получаем только те LD-значения, которые соответствуют
                                        #квадратикам, индекс строки которых больше индекса столбца.
                                        if row_index > col_index:
                                                
                                                #Обращение к оффлайн-калькулятору для
                                                #получения словаря с r2 и D' текущей пары SNP.
                                                #0-ячейка будет заменена на найденное значение
                                                #LD выбранной величины, либо останется нулевой.
                                                #Ноль сохранится, если пара находится в полном
                                                #равновесии по сцеплению, либо если LD этой
                                                #пары ниже заданного пользователем порога.
                                                #Значение LD округляется до 5 знаков.
                                                #Если это вам кажется неразумным, пишите в Issues.
                                                ld_vals = ld_calc(rows_by_chrs[chr_num][row_index],
                                                                  rows_by_chrs[chr_num][col_index],
                                                                  sample_indices)
                                                if 'ld_value_thres' in locals():
                                                        if ld_vals[thres_ld_measure] < ld_value_thres:
                                                                continue
                                                ld_two_dim[row_index][col_index] = round(ld_vals[ld_measure], 5)
                                                
                        ##Сохранение результатов.
                        
                        #Пользователь указал создавать текстовые версии LD-матриц.
                        if output == 'table' or output == 'both' or output == '':
                                
                                print(f'{src_file_name}: хромосома {chr_num}: сохранение текстовой LD-матрицы...')
                                
                                #Создание текстового конечного файла.
                                #Прописываем в него шапку, а потом строки,
                                #добавляя к каждой из них соответствующий refSNPID.
                                tsv_file_name = f'{src_file_base}_chr{chr_num}_LD_tab.tsv'
                                with open(os.path.join(trg_dir_path, tsv_file_name), 'w') as tsv_file_opened:
                                        tsv_file_opened.write('\t' + '\t'.join(rs_ids_srtd) + '\n')
                                        rs_id_index = 0
                                        for row in ld_two_dim:
                                                line = '\t'.join([str(cell) for cell in row]) + '\n'
                                                tsv_file_opened.write(rs_ids_srtd[rs_id_index] + '\t' + line)
                                                rs_id_index += 1
                                                
                        #Если переменная, ссылающаяся на значение, которое
                        #должно стать одним из аргументов функции построения
                        #диаграммы, существует, то это означает, что
                        #пользователь предпочёл LD-матрицы визуализировать.
                        if 'color_map' in locals():
                                
                                print(f'{src_file_name}: хромосома {chr_num}: визуализация LD-матрицы...')
                                
                                #Создание объекта диаграммы с помощью Plotly.
                                #Он представляет собой структуру со свойствами
                                #словаря, в которую вложены структуры со свойствами
                                #словарей и списков, содержащие параметры тепловой карты.
                                #С помощью аргументов функции create_annotated_heatmap
                                #можно менять настройки диаграммы, находящиеся
                                #в значении ключа верхнего уровня 'data'.
                                #Аргумент этой функции, который имеет возможность
                                #выбрать пользователь - цветовая схема тепловой карты.
                                #Не изменяемые в рамках интерактивного диалога аргументы:
                                #массивы с LD-данными и лейблами по X, Y, наличие
                                #разделительных линий между квадратиками, а также более
                                #тёмная заливка квадратиков при больших значениях в них.
                                ld_heatmap = ff.create_annotated_heatmap(ld_two_dim,
                                                                         x=rs_ids_srtd,
                                                                         y=rs_ids_srtd,
                                                                         xgap=1,
                                                                         ygap=1,
                                                                         colorscale=color_map,
                                                                         reversescale=True)
                                
                                #В аннотированных тепловых картах Plotly
                                #ось Y по умолчанию направлена снизу вверх.
                                #Чтобы ориентация тепловых карт, создаваемых
                                #этой программой, была такая же, как в
                                #LDMatrix, переворачиваем диаграмму по Y.
                                #Параметры осей заложены в структуры,
                                #принадлежащие ключу верхнего уровня 'layout'.
                                ld_heatmap['layout']['yaxis']['autorange'] = 'reversed'
                                
                                #Эта программа может опционально размещать
                                #значения LD в квадратиках тепловой карты,
                                #а также наносить лейблы на оси X и Y.
                                #Допустим, пользователь предпочёл не выводить
                                #на диаграмму никаких текстовых данных.
                                #Тогда, чтобы не вписывать ничего в квадратики, нужно
                                #сделать пустым списком значение ключа 'annotations',
                                #содержащего все настройки надписей для квадратиков,
                                #а в настройки осей, находящиеся в глубинах
                                #'layout', следует добавить ключ 'showticklabels'
                                #со значением, запрещающим вывод лейблов.
                                if texts == 'no' or texts == 'n' or texts == '':
                                        ld_heatmap['layout']['annotations'] = []
                                        ld_heatmap['layout']['xaxis']['showticklabels'] = False
                                        ld_heatmap['layout']['yaxis']['showticklabels'] = False
                                        
                                #Пользователь дал добро выводить на диаграмму надписи.
                                else:
                                        
                                        #Пользователь кастомизировал размер шрифта значений внутри квадратиков.
                                        #Тогда в каждый подсловарь структуры 'annotations' будет
                                        #добавлена пара ключ-значение с данным размером.
                                        if val_font_size != 'default':
                                                for ann_num in range(len(ld_heatmap['layout']['annotations'])):
                                                        ld_heatmap['layout']['annotations'][ann_num]['font']['size'] = int(val_font_size)
                                                        
                                        #По желанию, пользователь может задать недефолтное
                                        #значение параметра размера шрифта лейблов осей.
                                        if lab_font_size != 'default':
                                                ld_heatmap['layout']['xaxis']['tickfont'] = {'size': int(lab_font_size)}
                                                ld_heatmap['layout']['yaxis']['tickfont'] = {'size': int(lab_font_size)}
                                        
                                #Построение диаграммы и её сохранение в HTML.
                                html_file_path = f'{os.path.join(trg_dir_path, src_file_base)}_chr{chr_num}_LD_diag.html'
                                py.offline.plot(ld_heatmap, filename = html_file_path, auto_open=False)
                                
#Если задать для первой тепловой карты одну цветовую
#схему, а для второй - другую, то наблюдается баг:
#вторая диаграмма рисуется тех же цветов, что и первая.
#Лечится ежеразовым удалением питоновского кэша,
#образующегося при каждом запуске программы.
shutil.rmtree(os.path.join(os.getcwd(), 'backend', '__pycache__'))
