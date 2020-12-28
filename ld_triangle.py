__version__ = 'V11.3'

class PrepSingleProc():
        '''
        Класс, спроектированный
        под безопасное параллельное
        построение матриц значений
        неравновесия по сцеплению.
        '''
        def __init__(self, args):
                '''
                Получение атрибутов, необходимых заточенной
                под многопроцессовое выполнение функции
                формирования графических и/или текстовых
                LD-матриц. Атрибуты должны быть созданы
                единожды и далее ни в коем случае не
                изменяться. Получаются они в основном
                из указанных исследователем опций.
                '''
                self.src_dir_path = os.path.normpath(args.src_dir_path)
                self.intgen_dir_path = os.path.normpath(args.intgen_dir_path)
                if args.trg_top_dir_path is None:
                        self.trg_top_dir_path = self.src_dir_path
                else:
                        self.trg_top_dir_path = os.path.normpath(args.trg_top_dir_path)
                self.meta_lines_quan = args.meta_lines_quan
                if args.skip_intgen_data_ver:
                        self.intgen_convdb_path = os.path.join(self.intgen_dir_path,
                                                               'conversion.db')
                else:
                        self.intgen_convdb_path = prep_intgen_data(self.intgen_dir_path)
                if args.gend_names == 'male':
                        self.gend_names = ('male',)
                elif args.gend_names == 'female':
                        self.gend_names = ('female',)
                else:
                        self.gend_names = ('male', 'female')
                self.pop_names = tuple(args.pop_names.upper().split(','))
                self.sample_names = get_sample_names(self.gend_names,
                                                     self.pop_names,
                                                     self.intgen_convdb_path)
                self.ld_measure = args.ld_measure
                self.ld_low_thres = args.ld_low_thres
                self.matrix_type = args.matrix_type
                self.heatmap_json = args.heatmap_json
                self.disp_letters = args.disp_letters
                self.color_pal = args.color_pal
                self.font_size = args.font_size
                self.square_shape = args.square_shape
                self.dont_disp_footer = args.dont_disp_footer
                
        def create_matrix(self, src_file_name):
                '''
                Функция создания одной LD-матрицы.
                '''
                
                #Считывание исходной таблицы, извлечение оттуда
                #rsIDs и создание словаря, в котором позиции и
                #идентификаторы вариантов разбиты по хромосомам.
                data_by_chrs = create_src_dict(self.src_dir_path,
                                               src_file_name,
                                               self.meta_lines_quan,
                                               self.intgen_convdb_path)
                
                #В одну папку второго уровня планируется размещать все
                #результаты, полученные по данным одного исходного файла.
                src_file_base = src_file_name.rsplit('.', maxsplit=1)[0]
                trg_dir_path = os.path.join(self.trg_top_dir_path,
                                            f'{src_file_base}_LD_matr')
                
                #Для вариантов одной хромосомы
                #создастся одна матрица.
                for chrom in data_by_chrs:
                        
                        #Проверяем, набралось ли хотя
                        #бы 2 варианта, относящиеся к
                        #текущей хромосоме. Если да, то
                        #появляется смысл в реальном создании
                        #конечной папки второго уровня.
                        if len(data_by_chrs[chrom]) < 2:
                                continue
                        if os.path.exists(trg_dir_path) == False:
                                os.mkdir(trg_dir_path)
                                
                        #Чтобы потом проще было визуально оценивать
                        #влияние физического расстояния на LD,
                        #rsIDs отсортируются по геномным позициям.
                        data_by_chrs[chrom].sort(key=lambda row: row[0])
                        poss_srtd, rs_ids_srtd = [], []
                        for row in data_by_chrs[chrom]:
                                poss_srtd.append(row[0])
                                rs_ids_srtd.append(row[1])
                                
                        #Знание количества rsIDs в
                        #ближайшей перспективе пригодится,
                        #чтобы задать размеры матрицы,
                        #а в дальнейшей - чтобы
                        #оформить её табличную версию.
                        vars_quan = len(rs_ids_srtd)
                        
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
                        
                        #Построение шаблона квадратного двумерного массива, состоящего
                        #из нулей. Нули в дальнейшем могут заменяться на значения LD.
                        ld_two_dim = [[0 for col_index in range(vars_quan)] for row_index in range(vars_quan)]
                        
                        #В случае, если будет рисоваться диаграмма,
                        #такой же шаблон понадобится для создания
                        #матрицы сопутствующей информации.
                        if self.matrix_type in ['heatmap', 'both']:
                                info_two_dim = copy.deepcopy(ld_two_dim)
                                
                        #Для расчёта LD и аннотирования вариантов
                        #потребуются данные проекта 1000 Genomes.
                        #Открываем соответствующий текущей хромосоме
                        #tabix-индексированный 1000 Genomes-архив с
                        #помощью pysam. Pysam тут пригождается для
                        #быстрого доступа к случайным строкам архива.
                        with VariantFile(os.path.join(self.intgen_dir_path,
                                                      f'{chrom}.vcf.gz')) as intgen_vcf_opened:
                                
                                #Перебор индексов строк и столбцов
                                #изначально нулевых матриц.
                                for row_index in range(vars_quan):
                                        for col_index in range(vars_quan):
                                                
                                                #Матрица, в принципе, может
                                                #быть квадратом, состоящим
                                                #из двух одинаковых по форме
                                                #и содержимому прямоугольных
                                                #треугольников, разделённых
                                                #диагональю 0-ячеек. Думаю,
                                                #разумнее оставить лишь один
                                                #из этих треугольников. Для
                                                #этого получаем только те
                                                #значения, которые соответствуют
                                                #ячейкам двумерного массива,
                                                #индекс строки которых
                                                #больше индекса столбца.
                                                if row_index <= col_index:
                                                        continue
                                                
                                                #Вытаскивание из 1000 Genomes и отбор
                                                #по сэмплам фазированных генотипов текущей
                                                #пары вариантов. Разбиение пар генотипов
                                                #на отдельные, что необходимо из-за требования
                                                #калькулятора LD. Извлечение из 1000 Genomes
                                                #аннотаций каждого варианта обрабатываемой пары.
                                                y_var_genotypes, x_var_genotypes = [], []
                                                y_var_row = data_by_chrs[chrom][row_index]
                                                for intgen_rec in intgen_vcf_opened.fetch(chrom,
                                                                                          y_var_row[0] - 1,
                                                                                          y_var_row[0]):
                                                        if intgen_rec.id != y_var_row[1]:
                                                                continue
                                                        y_var_alleles = intgen_rec.ref + '/' + intgen_rec.alts[0]
                                                        y_var_type = intgen_rec.info['VT'][0]
                                                        for sample_name in self.sample_names:
                                                                try:
                                                                        y_var_genotypes += intgen_rec.samples[sample_name]['GT']
                                                                except KeyError:
                                                                        continue
                                                        break
                                                x_var_row = data_by_chrs[chrom][col_index]
                                                for intgen_rec in intgen_vcf_opened.fetch(chrom,
                                                                                          x_var_row[0] - 1,
                                                                                          x_var_row[0]):
                                                        if intgen_rec.id != x_var_row[1]:
                                                                continue
                                                        x_var_alleles = intgen_rec.ref + '/' + intgen_rec.alts[0]
                                                        x_var_type = intgen_rec.info['VT'][0]
                                                        for sample_name in self.sample_names:
                                                                try:
                                                                        x_var_genotypes += intgen_rec.samples[sample_name]['GT']
                                                                except KeyError:
                                                                        continue
                                                        break
                                                
                                                #Обращение к оффлайн-калькулятору
                                                #для получения словаря с r2, D' и
                                                #частотами альтернативных аллелей
                                                #пары вариантов для выбранных
                                                #исследователем популяций и полов.
                                                trg_vals = calc_ld(y_var_genotypes,
                                                                   x_var_genotypes)
                                                
                                                #Каждый элемент визуализируемой матрицы
                                                #аннотируется: параллельно с накоплением
                                                #массива LD-значений растёт массив дополнительной
                                                #информации по каждой паре вариантов.
                                                if self.matrix_type in ['heatmap', 'both']:
                                                        info_two_dim[row_index][col_index] = f'''
r2: {trg_vals["r_square"]}<br>
D': {trg_vals["d_prime"]}<br>
abs_dist: {abs(poss_srtd[col_index] - poss_srtd[row_index])}<br><br>
{rs_ids_srtd[col_index]}.hg38_pos: {poss_srtd[col_index]}<br>
{rs_ids_srtd[row_index]}.hg38_pos: {poss_srtd[row_index]}<br><br>
{rs_ids_srtd[col_index]}.alleles: {x_var_alleles}<br>
{rs_ids_srtd[row_index]}.alleles: {y_var_alleles}<br><br>
{rs_ids_srtd[col_index]}.type: {x_var_type}<br>
{rs_ids_srtd[row_index]}.type: {y_var_type}<br><br>
{rs_ids_srtd[col_index]}.alt_freq: {trg_vals['var_2_alt_freq']}<br>
{rs_ids_srtd[row_index]}.alt_freq: {trg_vals['var_1_alt_freq']}
'''
                                                        
                                                #Исследователь мог установить нижний порог LD.
                                                #Соответствующий блок кода неспроста расположен после
                                                #блока накопления аннотаций: на диаграммах клеточки с
                                                #подпороговыми LD будут закрашены как нулевые, но
                                                #зато при наведении курсора там отобразятся настоящие
                                                #LD-значения, как раз извлекаемые из массива с аннотациями.
                                                #При обратном расположении этих блоков аннотации подпороговых
                                                #LD не сохранялись бы, ведь в блоке фильтрации - continue.
                                                if self.ld_low_thres != None:
                                                        if trg_vals[self.ld_measure] < self.ld_low_thres:
                                                                continue
                                                        
                                                #Если значение LD не отсеилось как подпороговое,
                                                #то попадёт в LD-матрицу: 0-ячейка будет заменена
                                                #на найденное значение LD выбранной величины.
                                                ld_two_dim[row_index][col_index] = trg_vals[self.ld_measure]
                                                
                        #Стремящееся быть информативным название
                        #конечного файла. Какое к нему далее будет
                        #пристыковано расширение - зависит от
                        #выбранного исследователем формата.
                        trg_file_base = f'{src_file_base}_chr{chrom}_{self.ld_measure[0]}'
                        
                        #Визуализация матрицы с помощью plotly.
                        if self.matrix_type in ['heatmap', 'both']:
                                
                                #Исследователь дал добро
                                #выводить на диаграмму надписи:
                                #rsIDs в качестве лейблов осей и
                                #значения LD внутри квадратиков
                                #непосредственно тепловой карты.
                                if self.disp_letters:
                                        
                                        #Создание объекта аннотированной тепловой карты.
                                        #Из чего он состоит - см. в ридми к ld_triangle.
                                        #Здесь только отмечу, что create_annotated_heatmap -
                                        #высокоуровневая функция библиотеки plotly,
                                        #берущая на себя большую часть этой работы.
                                        ld_heatmap = ff.create_annotated_heatmap(ld_two_dim,
                                                                                 x=rs_ids_srtd,
                                                                                 y=rs_ids_srtd,
                                                                                 hovertext=info_two_dim,
                                                                                 hoverinfo='text',
                                                                                 xgap=1,
                                                                                 ygap=1,
                                                                                 colorscale=self.color_pal,
                                                                                 showscale=False)
                                        
                                        #Возможная кастомизация размера шрифта
                                        #подписей к осям и чисел в квадратиках.
                                        if self.font_size != None:
                                                ld_heatmap.layout.xaxis.tickfont.size = self.font_size
                                                ld_heatmap.layout.yaxis.tickfont.size = self.font_size
                                                for ann_num in range(len(ld_heatmap.layout.annotations)):
                                                        ld_heatmap['layout']['annotations'][ann_num]['font']['size'] = self.font_size
                                                        
                                #Исследователь предпочёл выводить на тепловую
                                #карту минимум текстовых данных. Жертвовать
                                #надписями обычно приходится во избежание их
                                #взамного наползания в крупных диаграммах.
                                #Построим объект смысловой части диаграммы
                                #и объект вторичных настроек, собирём
                                #их в финальный объект. Подробнее о
                                #структуре объектов plotly - в ридми.
                                else:
                                        trace = go.Heatmap(z=ld_two_dim,
                                                           hovertext=info_two_dim,
                                                           hoverinfo='text',
                                                           xgap=1,
                                                           ygap=1,
                                                           colorscale=self.color_pal,
                                                           showscale=False)
                                        layout = go.Layout(xaxis_showticklabels=False,
                                                           yaxis_showticklabels=False)
                                        ld_heatmap = go.Figure(data=trace,
                                                               layout=layout)
                                        
                                #Опциональное приведение
                                #диаграммы к квадратной форме.
                                if self.square_shape:
                                        ld_heatmap.update_layout(xaxis_constraintoward='left',
                                                                 yaxis_scaleanchor='x',
                                                                 yaxis_scaleratio = 1,
                                                                 plot_bgcolor='rgba(0,0,0,0)')
                                        
                                #Следующие настройки будут
                                #касаться, в основном, надписей,
                                #отличных от LD-значений и
                                #rsIDs - заголовка и футера.
                                #Чтобы размещать футер, пришлось
                                #пойти на небольшую хитрость -
                                #выделить под него тайтл оси X.
                                #Помимо всего прочего, переворачиваем
                                #диаграмму по Y ради визуальной
                                #совместимости с хитмэпами LDmatrix.
                                title = f'''
defines color: {self.ld_measure} ░
LD threshold: {self.ld_low_thres} ░
chromosome: {chrom} ░
genders: {", ".join(self.gend_names)} ░
populations: {", ".join(self.pop_names)}
'''
                                ld_heatmap.update_layout(title_text=title,
                                                         xaxis_side='bottom',
                                                         yaxis_autorange='reversed')
                                if self.dont_disp_footer == False:
                                        footer = '''
made by ld_triangle from <a href="https://github.com/PlatonB/ld-tools">ld-tools</a> ░
readme:
<a href="https://github.com/PlatonB/ld-tools/blob/master/README.md">ru</a>
<a href="https://github.com/PlatonB/ld-tools/blob/master/README-EN.md">en</a> ░
<a href="https://www.tinkoff.ru/rm/bykadorov.platon1/7tX2Y99140/">donate</a>
'''
                                        ld_heatmap.update_layout(xaxis_title_text=footer,
                                                                 xaxis_title_font_size=10)
                                        
                                #Прописывание всех данных диаграммы в
                                #JSON, если это необходимо исследователю.
                                if self.heatmap_json:
                                        debug_file_name = trg_file_base + '.json'
                                        ld_heatmap.write_json(os.path.join(trg_dir_path, debug_file_name),
                                                              pretty=True)
                                        
                                #Сохранение диаграммы в HTML.
                                html_file_name = trg_file_base + '.html'
                                ld_heatmap.write_html(os.path.join(trg_dir_path, html_file_name))
                                
                        #Исследователь выбрал опцию создавать
                        #табличные варианты LD-матриц.
                        if self.matrix_type in ['table', 'both']:
                                
                                #Создание текстового конечного файла. Прописываем в него хэдер
                                #с общими характеристиками матрицы, пустую строку и две шапки:
                                #одна - с rsIDs, другая - с позициями. Потом прописываем
                                #LD-строки, добавляя перед каждой из них тоже rsID и позицию.
                                tsv_file_name = trg_file_base + '.tsv'
                                with open(os.path.join(trg_dir_path, tsv_file_name), 'w') as tsv_file_opened:
                                        tab, poss_srtd = '\t', list(map(lambda pos: str(pos), poss_srtd))
                                        tsv_file_opened.write(f'##General\tinfo:\t{self.ld_measure}\tchr{chrom}\t{tab.join(self.pop_names)}\t{tab.join(self.gend_names)}\n\n')
                                        tsv_file_opened.write('rsIDs\t\t' + '\t'.join(rs_ids_srtd) + '\n')
                                        tsv_file_opened.write('\tPositions\t' + '\t'.join(poss_srtd) + '\n')
                                        for row_index in range(vars_quan):
                                                line = '\t'.join(map(str, ld_two_dim[row_index])) + '\n'
                                                tsv_file_opened.write(rs_ids_srtd[row_index] + '\t' +
                                                                      poss_srtd[row_index] + '\t' +
                                                                      line)
                                                
####################################################################################################

import sys, locale, os, datetime, copy

#Подавление формирования питоновского кэша с
#целью предотвращения искажения результатов.
sys.dont_write_bytecode = True

from cli.ld_triangle_cli_ru import add_args_ru
from cli.ld_triangle_cli_en import add_args_en
from backend.prep_intgen_data import prep_intgen_data
from backend.get_sample_names import get_sample_names
from backend.create_src_dict import create_src_dict
from multiprocessing import Pool
from pysam import VariantFile
from backend.calc_ld import calc_ld
import plotly.figure_factory as ff
import plotly.graph_objs as go

#Подготовительный этап: обработка
#аргументов командной строки,
#создание экземпляра содержащего
#ключевую функцию класса, определение
#оптимального количества процессов.
if locale.getdefaultlocale()[0][:2] == 'ru':
        args = add_args_ru(__version__)
else:
        args = add_args_en(__version__)
prep_single_proc = PrepSingleProc(args)
max_proc_quan = args.max_proc_quan
src_file_names = os.listdir(prep_single_proc.src_dir_path)
src_files_quan = len(src_file_names)
if max_proc_quan > src_files_quan <= 8:
        proc_quan = src_files_quan
elif max_proc_quan > 8:
        proc_quan = 8
else:
        proc_quan = max_proc_quan
        
print('\nLD matrix(-es) creation')
print(f'\tnumber of parallel processes: {proc_quan}')

#Параллельный запуск создания матриц. Замер времени
#выполнения этого кода с точностью до микросекунды.
with Pool(proc_quan) as pool_obj:
        exec_time_start = datetime.datetime.now()
        pool_obj.map(prep_single_proc.create_matrix, src_file_names)
        exec_time = datetime.datetime.now() - exec_time_start
        
print(f'\tparallel computation time: {exec_time}')
