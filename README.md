# Синопсис.

## Подготовка.
1. Установите [IDLE](https://github.com/PlatonB/bioinformatic-python-scripts#установка-среды-разработки) или другую удобную вам среду разработки, поддерживающую Python 3. Если у вас Linux, можете, не используя среду разработки, запускать Python-программы из терминала.
2. Установите необходимые [Python-модули](https://github.com/PlatonB/ld-tools#установка-сторонних-python-модулей).
3. Скачайте архив с программой:
```
Clone or download
```
(зелёная кнопка наверху страницы репозитория).

4. Распакуйте его в любую папку.

### Установка сторонних Python-модулей.
#### [Linux](https://elementary.io/ru/).
```
sudo pip install название_модуля
```

#### Windows.
TBD

#### Модули.
- Для компонента, отвечающего за поиск SNPs в LD с запрашиваемыми, понадобится `pysam` — Python-обёртка к Samtools.
- Для компонента визуализации LD-матриц нужен будет `plotly`.

## Запуск.
1. Если пользуетесь `IDLE`, откройте в нём нужный компонент программы и запустите его кнопкой `F5`.
2. Введите [опции](https://github.com/PlatonB/bioinformatic-python-scripts#диалог-с-пользователем) в рамках интерактивного диалога.

В `IDLE`, когда работа программы завершится, появятся символы `>>>`.

# LD_triangle.
Программа LD_triangle строит матрицы значений неравновесия по сцеплению для пар SNP. Каждая матрица создаётся по всем парным сочетаниям соответствующего единого набора refSNPID. Важно отметить, что если уже найдено LD для пары SNPx-SNPy, то обработка пары SNPy-SNPx допущена не будет. Поэтому вместо квадратов с результатами получаются прямоугольные треугольники. Программа является написанной с нуля десктопной реализацией идеи компонента LDmatrix веб-сервиса LDlink.

## Ключевые отличия от [LDmatrix](https://ldlink.nci.nih.gov/?tab=ldmatrix) 3.4.0:
- Работа с исходными файлами:
	- не нужно вручную загружать по одному файлу, программа сама сканирует всю папку с данными
	- автоматически распознаётся столбец с идентификаторами SNP каждой таблицы;
- После подготовки кэша работу можно проводить полностью оффлайн;
- Возможность тонкой настройки с помощью дружественного интерактивного диалога;
- Не возникает конфликт, если на вход подаются SNP разных хромосом: для каждой хромосомы создастся отдельная LD-матрица;
- Возможность обнулить те значения LD, которые ниже определённого порога, чтобы при рассмотрении матрицы было проще сосредоточиться на наиболее сильно сцепленных парах;
- (ИМХО!) Более наглядные и симпатичные диаграммы (можете оспорить это в [Issues](https://github.com/PlatonB/ld-tools/issues):smile:)

## Результаты.
На пересечении строки с одним SNP и столбца с другим отображается значение LD. Оно может быть представлено в виде вещественного числа, степени насыщенности окраски квадратика, либо сочетанием одного и другого.
![Пример тепловой карты с указанием закономерностей физического расстояния](https://raw.githubusercontent.com/PlatonB/ld-tools/master/LD_triangle_output_with_explanations.png)
Программа сортирует SNP по возрастанию их позиции в хромосоме. В LD-матрице обозначающие SNP лейблы располагаются таким образом, что по Y увеличение позиции SNPs идёт сверху вниз, а по X - слева направо. Исходя из этого, самое маленькое физическое расстояние будет у SNP-пар на гипотенузе, а самое большое - около пересечения катетов. Если конкретнее, то относительно любого квадратика физическое расстояние внутри каждой пары увеличивается влево, вниз и по "результирующей" диагонали (см. рисунок).

## Корректировка диаграммы (для опытных пользователей).
```
Не закрывайте и не перезапускайте интерактивную консоль, пока не просмотрели полученную диаграмму!
```

Если после появления `>>>` консоль IDLE всё ещё открыта и не перезапущена вами (`CTRL+F6`), созданные в процессе работы программы объекты хранятся в оперативной памяти. Этим можно воспользоваться, пересоздав диаграмму без повторного выполнения целого скрипта. Если же запускать скрипт по новой, то придётся ждать, когда все значения LD заново сосчитаются, что может длиться очень долго.

Если вы знаете основы Python, то без труда разберётесь в объекте, формирующим диаграмму. Выведите этот объект на экран:
```
ld_heatmap
```

(дважды кликните на `squeezed text`).

К элементам этой структуры можно обращаться как к элементам обычных словарей и списков. Можно добавлять новые элементы и заменять имеющиеся. Разве что `del` не работает, но эту проблему можно обойти. Если вы - не Python-программист и не разбираетесь в структурах данных, то просто копируйте в интерактивную консоль приведённые ниже строки кода.

Самая частая неприятность - когда надписи наползают друг на друга. Если диаграмма не очень большая (в пределах 50х50), то проблема решается уменьшением шрифта. В Plotly размер шрифта по умолчанию - 12, поэтому в качестве примера уменьшенного шрифта здесь и далее привожу 10. Смело пробуйте и другие значения. Если что, объект диаграммы можно будет пересоздать (см. ниже), тем самым сбросив ваши изменения.

В первую очередь, по мере увеличения тепловой карты начинают пересекаться LD-значения внутри квадратиков. Сделаем их поменьше:
```
for ann_num in range(len(ld_heatmap['layout']['annotations'])):
	ld_heatmap['layout']['annotations'][ann_num]['font']['size'] = 10
```

После каждого внесённого изменения предпросматривайте тепловую карту (при этом будет сохраняться только временный файл):
```
py.offline.plot(ld_heatmap)
```

Если размер матрицы ещё больше, то могут сливаться лейблы осей. Компактизируем их:
```
ld_heatmap['layout']['xaxis']['tickfont'] = {'size': 10}
```
```
ld_heatmap['layout']['yaxis']['tickfont'] = {'size': 10}
```

Кое-как ещё умещаются лейблы осей, но совсем не влезают надписи в квадратиках? Тогда последние удалим. Как уже говорилось ранее, с помощью `del` это не произвести, но можно схитрить - вместо ненужного объекта создать пустой объект того же типа. Взяв на вооружение этот способ, упраздняем список с настройками текстового заполнения квадратиков, тем самым освобождая квадратики от надписей:
```
ld_heatmap['layout']['annotations'] = []
```

Если матрица настолько крупная, что лейблы осей тоже стали выходить из берегов, то удаляем и их:
```
ld_heatmap['layout']['xaxis']['showticklabels'] = False
```
```
ld_heatmap['layout']['yaxis']['showticklabels'] = False
```

Получится тепловая карта такого вида:
![Пример крупной тепловой карты без надписей](https://raw.githubusercontent.com/PlatonB/ld-tools/master/LD_triangle_big_heatmap.png)

Вновь что-то кажется кривым? Продолжайте описанные выше эксперименты с параметрами. Если же проделанная работа устраивает, можете сохранить обновлённую диаграмму в ту папку, путь к которой указывали в рамках пользовательского диалога. Учтите, что старая диаграмма при этом заменится на новую.
```
py.offline.plot(ld_heatmap, filename = html_file_path, auto_open=False)
```

Может так сложиться, что ваши эксперименты зашли слишком далеко и нужно вернуть тепловой карте тот вид, который она имела непосредственно после завершения работы программы. Тогда создадим объект диаграммы заново:
```
ld_heatmap = ff.create_annotated_heatmap(snp_pair_ld_vals, x=rs_ids_srtd, y=rs_ids_srtd, xgap=1, ygap=1, colorscale=color_map, reversescale=True)
```

Пересоздавая объект тепловой карты, можете заодно поиграть с аргументами соответствующей функции. К примеру, настройте ширину промежутков между квадратиками (`xgap=2`, `ygap=2`), поменяйте цветовую схему (`colorscale='Hot'`).

Перевернём тепловую карту по оси Y, чтобы она выглядела, как результат работы LDmatrix:
```
ld_heatmap['layout']['yaxis']['autorange'] = 'reversed'
```

Ещё для поддержания схожести с диаграммами LDmatrix можно сделать наш вариант полностью квадратным. Т.е. приведём к квадратной форме и элементы тепловой карты, и область построения (допустим, 1000х1000 пикселей).
```
ld_heatmap['layout']['yaxis']['scaleanchor'] = 'x'
```
```
ld_heatmap['layout']['yaxis']['scaleratio'] = 1
```
```
ld_heatmap['layout']['width'] = 1000
```
```
ld_heatmap['layout']['height'] = 1000
```

В этом примере созданная по новой диаграмма будет выглядеть как-то так:
![Пример тепловой карты 1000х1000 пикселей с квадратными элементами](https://raw.githubusercontent.com/PlatonB/ld-tools/master/LD_triangle_quadratic_heatmap.png)

Я охватил только самые напрашивающиеся идеи редактирования получаемых программой тепловых карт. Если вы обладаете колоссальным терпением и непреодолимым желанием производить тонкую настройку диаграмм Plotly, то эта ссылка для вас:
https://plot.ly/python/reference/
