# Синопсис.
## Компоненты.
| Программа | Предназначение |
| --------- | -------------- |
| ld_area | выводит для каждого из запрашиваемых SNPs все снипы в пределах двух порогов: не ниже определённого LD и не дальше заданных фланков |
| ld_triangle | строит матрицы значений LD всех возможных пар SNPs исходного набора |

## Преимущества.
До появления _ld_tools_ уже существовали мощные, популярные, разрабатываемые крупными лабораториями решения: _LDlink_, _Ensembl REST API_, _Haploreg_ и др.. На сегодняшний день по ряду функциональных и юзабилити-моментов _ld_tools_ занимает более выигрышную позицию.
- Не нужно вручную загружать по одному файлу, программа сама сканирует всю исходную папку;
- Автоматически распознаётся столбец с идентификаторами SNP каждой таблицы.
- Вплоть до 8 исходных файлов могут обрабатываться параллельно.
- После подготовки кэша работу можно проводить полностью оффлайн. Никаких тормознутых API и 503-сюрпризов.
- Распределение результатов по хромосомным подпапкам.
- Возможность тонкой настройки каждого запуска с помощью агрументов [командной строки](https://github.com/PlatonB/ngs-pipelines#преодолеваем-страх-командной-строки-linux).
- Сохранение в конечные файлы настроек их получения.
- Геномные координаты - только GRCh38. Нет характерной для биоинформатических тулзов путаницы по этому вопросу.
- Человеко-понятные имена переменных и подробные комментарии к коду. Делаю всё, чтобы проект пригождался в том числе и для изучения программирования в биоинформатике.
- Отличная производительность на неигровых ноутбуках. Тестировал на ASUS VivoBook X712FA с довольно хилым Intel Core i5-8265U.
- [Техподдержка](https://github.com/PlatonB/ld-tools/issues) лично от создателя программы, без посылания на Stack Overflow.

## Перед началом работы.
1. Хотя бы чуть-чуть ознакомьтесь с [командной строкой](https://github.com/PlatonB/ngs-pipelines#преодолеваем-страх-командной-строки-linux).
2. Установите сторонние Python-компоненты: необходимый для парсинга VCF-файлов _pysam_ и визуализатор _plotly_. Советую это делать посредством [_Conda_](https://github.com/PlatonB/ngs-pipelines#установка-conda).
```
conda install pysam=0.15.4 plotly
```

3. Скачайте архив с программой, нажав на зелёную кнопку наверху страницы репозитория.
```
Code
```
```
Download ZIP
```

4. Распакуйте его в любое место.
5. В терминале перейдите в разархивированную папку.
```
cd /path/to/ld-tools-master
```

6. Выведите справку по интересующему компоненту:
```
python3 ld_triangle -h
```

## Ограничения.
- Мультиаллельные SNPs, а также SNPs с не-rs-идентификаторами не поддерживаются, поэтому будут проигнорированы любым компонентом этой программы.
- Лучше исключать из исходных выборок данные X и Y-хромосом по причине нерешённой разработчиками 1000 Genomes [проблемы](https://github.com/samtools/bcftools/issues/1154) с PAR-регионами.
- Нельзя отключать интернет при подготовке данных, необходимых для быстрого оффлайн-вычисления LD. Этот процесс производится лишь однократно и занимает примерно полдня.
- Большие тепловые карты LD-значений, вроде 500х500 rsIDs, могут попросту не прорисоваться вашим браузером.

# LD_area.
Каждый исходный SNP порождает отдельный файл надпорогово-сцепленных с ним SNPs. Идея программы навеяна ныне заброшенным онлайн-инструментом _HaploReg_.

## Результаты.
- TSV. Формат по умолчанию. Это - обычные таблицы, удобные и информативные для клиницистов, в то же время, легко обрабатываемые [скриптами](https://github.com/PlatonB/bioinformatic-python-scripts). Ключевые столбцы - r2, D' и дистанция между найденными и запрашиваемой мутациями. Полезное дополнение - базовые аннотации каждой мутации (в том числе запрашиваемой): позиция, идентификатор, референсный и альтернативный аллели, тип, а также частота альтернативного аллеля.
![Пример TSV](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_area_TSV_output.png)
- JSON. По содержанию совпадает с TSV, но больше ориентирован на программистов.
![Пример JSON](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_area_JSON_output.png)
- Набор одних лишь rsIDs. Такая разновидность вывода лучше всего подходит для дальнейшего аннотирования мутаций специализированными программами ([_high_perf_bio_](https://github.com/PlatonB/high-perf-bio), _Ensembl VEP_ и др.).

# LD_triangle.
Каждая матрица создаётся по всем парным сочетаниям соответствующего единого набора rsIDs. Важно отметить, что если уже найдено LD для пары SNPx-SNPy, то обработка пары SNPy-SNPx допущена не будет. Поэтому вместо квадратов с результатами получаются прямоугольные треугольники. Программа является написанной с нуля десктопной реализацией идеи компонента _LDmatrix_ веб-сервиса _LDlink_.

## Ключевые особенности.
- Вывод в виде диаграмм и/или таблиц.
![Пример диаграммы](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_triangle_heatmap_example.png)
![Пример таблицы](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_triangle_text_example.png)

- Не возникает конфликт, если на вход подаются SNP разных хромосом: для каждой хромосомы создастся отдельная LD-матрица.
- Возможность обнулить те значения LD, которые ниже определённого порога, чтобы при рассмотрении матрицы было проще сосредоточиться на наиболее сильно сцепленных парах.
![Пример диаграммы без фильтрации](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_triangle_full_heatmap.png)
![Пример диаграммы с r2 >= 0.9](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_triangle_filtered_heatmap.png)

- Высокоинформативные всплывающие лейблы тепловых карт.
![Пример ховертекста](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_triangle_heatmap_with_hovertext.png)

- (ИМХО!) Более наглядные и симпатичные диаграммы, чем у LDmatrix (можете оспорить это в [Issues](https://github.com/PlatonB/ld-tools/issues):smile:).

## Подходящие цветовые палитры.
Подобрал такие палитры _Plotly_, с которыми LD-карты не будут смотреться вырвиглазно.

_algae, amp, blues, blugrn, bluyl, brwnyl, bugn, bupu, burg, burgyl, darkmint, deep, dense, emrld, gnbu, greens, greys, magenta, matter, mint, oranges, orrd, oryel, peach, pinkyl, pubu, pubugn, purd, purp, purples, purpor, rdpu, redor, reds, speed, sunset, sunsetdark, teal, tealgrn, tempo, turbid, ylgn, ylgnbu, ylorbr, ylorrd_

## Результаты.
На пересечении строки с одним SNP и столбца с другим отображается значение LD. Оно может быть представлено в виде вещественного числа, степени насыщенности окраски квадратика, либо сочетанием одного и другого.
![Пример тепловой карты с указанием закономерностей физического расстояния](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_triangle_annotated_heatmap_with_explanations.png)
Программа сортирует SNP по возрастанию их позиции в хромосоме. В LD-матрице обозначающие SNP лейблы располагаются таким образом, что по Y увеличение позиции SNPs идёт сверху вниз, а по X - слева направо. Исходя из этого, самое маленькое физическое расстояние будет у SNP-пар на гипотенузе, а самое большое - около пересечения катетов. Если конкретнее, то относительно любого квадратика физическое расстояние внутри каждой пары увеличивается влево, вниз и по "результирующей" диагонали (см. рисунок выше).

## Plotly.
Делать небольшие полуоффтопные тьюториалы - это уже традиция для моих ридми:). _ld_triangle_ использует библиотеку визуализации _plotly_, поэтому считаю нужным немного разобрать именно её. Это будет не огромная лавина отрывков кода для бездумного копирования, что привычно видеть на отечественных IT-порталах, а попытка лаконично акцентировать внимание на общем принципе. Поняв его, вы станете готовить графики без мучительного гугления.

Самый главный тезис - это что любая диаграмма _plotly_ строится по объекту класса `Figure`. Он состоит из вложенных друг в друга объектов со свойствами обычных питоновских словарей и списков. Если у вас есть хотя бы базовые знания _Python_, вы можете легко читать и редактировать `Figure`-объект. К тому же, понимание его строения позволяет быстро ориентироваться в длиннющей простыне [официальной документации](https://plot.ly/python/reference/).

Сохранить целиком объект тепловой карты попарных LD-значений можно флагом `-j` или `--heatmap-json` программы _ld_triangle_. Приводимые в примерах пути заменяйте на свои.
```
cd $HOME/Биоинформатика/ld-tools-master
```
```
python ld_triangle.py -S $HOME/Биоинформатика/SNP_data -D $HOME/Биоинформатика/1000_Genomes_data -f -e eur -j
```

Почему я написал "целиком"? Как правило, когда вы создаёте объект диаграммы, вам достаточно поменять пару-тройку параметров, а по остальным довериться дефолтам. Выводимый с помощью `print` отрывок объекта состоит из формирующих диаграмму данных и в подавляющей части недефолтных конфигурационных элементов. Но в JSON, помимо данных, отправляется абсолютно весь набор настроек. А если мы вручную извлечём содержимое JSON и выведем на экран, то увидим среди настроечных элементов преимущественно те, которые в явном виде были модифицированы на этапе создания диаграммы (в нашем случае внесением таких модификаций занимается _ld_triangle_). Теперь ещё раз, но короче: JSON содержит 100% кирпичиков объекта, а принт показывает только данные и в основном подкрученные настройки.
```
python3
```
```
from plotly.io import read_json
```
```
ld_heatmap = read_json('/home/yourname/Биоинформатика/результаты/SNPs_LD_matr/chr6_r_SNPs.json')
```
```
print(ld_heatmap)
```

То, что объект на самом деле считался из JSON в полном объёме, легко доказывается построением на его основе тепловой карты. Она получится такой же, как и изготовленная до этого полноценной программой.
```
ld_heatmap.show()
```

Теперь детально разберём отредактированную _ld_triangle_ часть объекта. Некоторые попавшие в принт элементы, правда, не были явно затронуты программой. Для простоты я сделал тепловую карту всего лишь по трём SNPs. Вот так будет выглядеть фрагмент её подкапотной части:
```
Figure({
    'data': [{'colorscale': [[0.0, 'rgb(247,252,245)'], [0.125,
                             'rgb(229,245,224)'], [0.25, 'rgb(199,233,192)'],
                             [0.375, 'rgb(161,217,155)'], [0.5,
                             'rgb(116,196,118)'], [0.625, 'rgb(65,171,93)'], [0.75,
                             'rgb(35,139,69)'], [0.875, 'rgb(0,109,44)'], [1.0,
                             'rgb(0,68,27)']],
              'hoverinfo': 'text',
              'hovertext': [[0, 0, 0], ["\nr2: 0.0003<br>\nD':
                            0.0247<br>\nabs_dist: 1060331<br><br>\nrs1521.hg38_pos:
                            31382927<br>\nrs8084.hg38_pos:
                            32443258<br><br>\nrs1521.alleles:
                            C/T<br>\nrs8084.alleles: A/C<br><br>\nrs1521.alt_freq:
                            0.7376<br>\nrs8084.alt_freq: 0.5865\n", 0, 0], ["\nr2:
                            0.0027<br>\nD': 0.0668<br>\nabs_dist:
                            1060942<br><br>\nrs1521.hg38_pos:
                            31382927<br>\nrs7192.hg38_pos:
                            32443869<br><br>\nrs1521.alleles:
                            C/T<br>\nrs7192.alleles: T/G<br><br>\nrs1521.alt_freq:
                            0.7376<br>\nrs7192.alt_freq: 0.6332\n", "\nr2:
                            0.8216<br>\nD': 1.0<br>\nabs_dist:
                            611<br><br>\nrs8084.hg38_pos:
                            32443258<br>\nrs7192.hg38_pos:
                            32443869<br><br>\nrs8084.alleles:
                            A/C<br>\nrs7192.alleles: T/G<br><br>\nrs8084.alt_freq:
                            0.5865<br>\nrs7192.alt_freq: 0.6332\n", 0]],
              'reversescale': False,
              'showscale': False,
              'type': 'heatmap',
              'x': [rs1521, rs8084, rs7192],
              'xgap': 1,
              'y': [rs1521, rs8084, rs7192],
              'ygap': 1,
              'z': [[0, 0, 0], [0.0003, 0, 0], [0.0027, 0.8216, 0]]}],
    'layout': {'annotations': [{'font': {'color': '#000000'},
                                'showarrow': False,
                                'text': '0',
                                'x': 'rs1521',
                                'xref': 'x',
                                'y': 'rs1521',
                                'yref': 'y'},
                               {'font': {'color': '#000000'},
                                'showarrow': False,
                                'text': '0',
                                'x': 'rs8084',
                                'xref': 'x',
                                'y': 'rs1521',
                                'yref': 'y'},
                               {'font': {'color': '#000000'},
                                'showarrow': False,
                                'text': '0',
                                'x': 'rs7192',
                                'xref': 'x',
                                'y': 'rs1521',
                                'yref': 'y'},
                               {'font': {'color': '#000000'},
                                'showarrow': False,
                                'text': '0.0003',
                                'x': 'rs1521',
                                'xref': 'x',
                                'y': 'rs8084',
                                'yref': 'y'},
                               {'font': {'color': '#000000'},
                                'showarrow': False,
                                'text': '0',
                                'x': 'rs8084',
                                'xref': 'x',
                                'y': 'rs8084',
                                'yref': 'y'},
                               {'font': {'color': '#000000'},
                                'showarrow': False,
                                'text': '0',
                                'x': 'rs7192',
                                'xref': 'x',
                                'y': 'rs8084',
                                'yref': 'y'},
                               {'font': {'color': '#000000'},
                                'showarrow': False,
                                'text': '0.0027',
                                'x': 'rs1521',
                                'xref': 'x',
                                'y': 'rs7192',
                                'yref': 'y'},
                               {'font': {'color': '#FFFFFF'},
                                'showarrow': False,
                                'text': '0.8216',
                                'x': 'rs8084',
                                'xref': 'x',
                                'y': 'rs7192',
                                'yref': 'y'},
                               {'font': {'color': '#000000'},
                                'showarrow': False,
                                'text': '0',
                                'x': 'rs7192',
                                'xref': 'x',
                                'y': 'rs7192',
                                'yref': 'y'}],
               'template': '...',
               'title': {'text': ('\ndefines color: r_square ░\nLD ' ... 'le, female ░\npopulations: EUR\n')},
               'xaxis': {'dtick': 1,
                         'gridcolor': 'rgb(0, 0, 0)',
                         'side': 'bottom',
                         'ticks': '',
                         'title': {'font': {'size': 10},
                                   'text': ('\nсделано с помощью ld_triangle' ... '6">поддержать разработчика</a>')}},
               'yaxis': {'autorange': 'reversed', 'dtick': 1, 'ticks': '', 'ticksuffix': '  '}}
})
```

Наиболее значимые компоненты словареподобного объекта `Figure` - `data` и `layout`. Первый посвящён самим данным и тесно с ними связанным настройкам, второй позволяет колдовать над дизайном диаграммы. `ld_heatmap['data'][0]['z']` - собственно, посчитанные для всех сниповых пар значения LD. `['data'][0]['x']` и `['data'][0]['y']` - надписи к осям, а именно, rsIDs. `['data'][0]['hovertext']` ведёт к двумерному массиву отформатированных с помощью HTML всплывающих характеристик каждой пары. `['layout']['annotations']` - вписываемые в квадратики значения LD и соответствующие настройки. `['layout']['title']['text']` - заголовок всего хитмэпа, куда помещаются биологические характеристики последнего. `['layout']['xaxis']['title']['text']` - заголовок оси X, внезапно приспособленный мною под футер с полезными ссылками. `['layout']['yaxis']['autorange']` при значении `reversed` переворачивает диаграмму по оси Y, что делается в нашем случае для унификации с хитмэпами LDmatrix.

В предыдущем абзаце, как видите, я ссылался на различные атрибуты чисто питоновским способом - обращением к словарным ключам и индексам элементов списков. Изменять атрибуты можно таким же макаром. Давайте, к примеру, компактизируем лейблы осей.
```
ld_heatmap['layout']['xaxis']['tickfont']['size'] = 10
```
```
ld_heatmap['layout']['yaxis']['tickfont']['size'] = 10
```

Экземпляр класса _Figure_ разрешает нам также задействовать ООП-синтаксис, чем не могут похвастать обычные словари.
```
ld_heatmap.layout.xaxis.tickfont.size = 10
```
```
ld_heatmap.layout.yaxis.tickfont.size = 10
```

Есть ещё замечательный синтаксический сахар - функция обновления `layout`.
```
ld_heatmap.update_layout(xaxis_tickfont_size=10)
```
```
ld_heatmap.update_layout(yaxis_tickfont_size=10)
```

Будьте осторожны: редактирование более высокоуровневых слоёв приводит к сбросу более низкоуровневых до значений по умолчанию. Например, вот так не делайте, иначе потеряете все отличные от кегля лейблов твики осей:
```
ld_heatmap.layout.xaxis = {'tickfont': {'size': 10}}
```
```
ld_heatmap.layout.yaxis = {'tickfont': {'size': 10}}
```

Послесловие. Любые атрибуты можно точно также легко доставать и править. При внимательном разглядывании полной редакции того или иного plotly-объекта можно даже совсем забыть про документацию и форумные треды.
