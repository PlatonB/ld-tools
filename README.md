# Синопсис.
## Преимущества.
Основанное на собственном опыте сравнение LD_tools с такими известными решениями, как LDlink, Ensembl REST API и Haploreg.
- Работа с исходными файлами:
	- не нужно вручную загружать по одному файлу, программа сама сканирует всю папку с данными;
	- автоматически распознаётся столбец с идентификаторами SNP каждой таблицы.
- После подготовки кэша работу можно проводить полностью оффлайн. Никаких тормознутых API и 503-сюрпризов.
- Возможность тонкой настройки с помощью дружественного интерактивного диалога. Не требуется командная строка.
- Геномные координаты - только GRCh38. Нет характерной для биоинформатических тулзов путаницы по этому вопросу.
- Код на чистом Python с человеко-понятными именами переменных и подробными комментариями. Делаю всё, чтобы проект пригождался в том числе и для изучения программирования в биоинформатике.

## Перед началом работы.
1. Установите [IDLE](https://github.com/PlatonB/bioinformatic-python-scripts#установка-среды-разработки) или другую удобную вам среду разработки, поддерживающую Python 3. Если у вас Linux, можете, не используя среду разработки, запускать Python-программы из терминала.
2. Установите необходимые [Python-модули](https://github.com/PlatonB/ld-tools#установка-сторонних-python-модулей).
3. Скачайте архив с программой:
```
Clone or download
```
(зелёная кнопка наверху страницы репозитория).

4. Распакуйте его в любую папку.
5. Создайте папку для данных из 1000 Genomes.

### Установка сторонних Python-модулей.
- Для всех компонентов потребуется `plyvel` (Python-обёртка к LevelDB).
- Для LD_area — компонента, отвечающего за поиск SNPs в LD с запрашиваемыми — понадобится `pysam` (Python-обёртка к Samtools).
- Для LD_triangle — компонента построения LD-матриц — нужны будут `plotly` (библиотека визуализации) и `numpy` (математический пакет в качестве обязательной зависимости для plotly).

#### [Linux](https://elementary.io/ru/).
```
sudo pip3 install plyvel pysam plotly numpy
```

#### Windows.
В Windows пока не работает `pysam`, поэтому запускать удастся только LD_triangle.
```
pip3 install plyvel plotly numpy
```

## Запуск.
1. Если пользуетесь `IDLE`, откройте в нём нужный компонент программы и запустите его кнопкой `F5`.
2. Введите [опции](https://github.com/PlatonB/bioinformatic-python-scripts#диалог-с-пользователем) в рамках интерактивного диалога.

В `IDLE`, когда работа программы завершится, появятся символы `>>>`.

## Важно!
- Мультиаллельные SNP не поддерживаются, поэтому будут проигнорированы любым компонентом этой программы.
- Подготовка данных, необходимых для быстрого оффлайн-вычисления LD, может занимать длительное время и сопровождаться небольшими лагами ОС. Но выполняться она будет только при первом запуске программы. Не отключайте интернет, и просто перетерпите сутки-двое:).
- По любым вопросам пишите в [Issues](https://github.com/PlatonB/ld-tools/issues) (потребуется регистрация на Github).
- Другие полезные биоинформатические инструменты моей разработки можете найти в репозитории https://github.com/PlatonB/bioinformatic-python-scripts.

# LD_area.
Выводит для каждого из запрашиваемых SNP все SNP в пределах двух порогов: не ниже определённого LD и не дальше заданных фланков.

## Ключевые особенности.
- Файлы, получающиеся на выходе:
	- JSON, содержащие значения LD и физической дистанции, а также идентификаторы и наиболее важные характеристики каждого SNP.
![Пример JSON](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_area_JSON_output.png)
	- Наборы одних лишь refSNPID для дальнейшего аннотирования специализированными программами (SnpEff, Ensembl VEP и др.).
- Распределение результатов по хромосомным подпапкам.

# LD_triangle.
Строит матрицы значений неравновесия по сцеплению для пар SNP. Каждая матрица создаётся по всем парным сочетаниям соответствующего единого набора refSNPID. Важно отметить, что если уже найдено LD для пары SNPx-SNPy, то обработка пары SNPy-SNPx допущена не будет. Поэтому вместо квадратов с результатами получаются прямоугольные треугольники. Программа является написанной с нуля десктопной реализацией идеи компонента LDmatrix веб-сервиса LDlink.

## Ключевые особенности.
- Вывод в виде диаграмм и/или таблиц.
![Пример диаграммы](https://raw.githubusercontent.com/PlatonB/ld-tools/master/LD_triangle_heatmap_output.png)
![Пример таблицы](https://raw.githubusercontent.com/PlatonB/ld-tools/master/LD_triangle_text_output.png)

- Не возникает конфликт, если на вход подаются SNP разных хромосом: для каждой хромосомы создастся отдельная LD-матрица.
- Возможность обнулить те значения LD, которые ниже определённого порога, чтобы при рассмотрении матрицы было проще сосредоточиться на наиболее сильно сцепленных парах.
![Пример диаграммы без фильтрации](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_triangle_full_heatmap.png)
![Пример диаграммы с r2 >= 0.9](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_triangle_filtered_heatmap.png)

- Высокоинформативные всплывающие лейблы тепловых карт.
![Пример ховертекста](https://raw.githubusercontent.com/PlatonB/ld-tools/master/LD_triangle_heatmap_hovertext.png)

- (ИМХО!) Более наглядные и симпатичные диаграммы (можете оспорить это в [Issues](https://github.com/PlatonB/ld-tools/issues):smile:).

## Результаты.
На пересечении строки с одним SNP и столбца с другим отображается значение LD. Оно может быть представлено в виде вещественного числа, степени насыщенности окраски квадратика, либо сочетанием одного и другого.
![Пример тепловой карты с указанием закономерностей физического расстояния](https://raw.githubusercontent.com/PlatonB/ld-tools/master/LD_triangle_output_with_explanations.png)
Программа сортирует SNP по возрастанию их позиции в хромосоме. В LD-матрице обозначающие SNP лейблы располагаются таким образом, что по Y увеличение позиции SNPs идёт сверху вниз, а по X - слева направо. Исходя из этого, самое маленькое физическое расстояние будет у SNP-пар на гипотенузе, а самое большое - около пересечения катетов. Если конкретнее, то относительно любого квадратика физическое расстояние внутри каждой пары увеличивается влево, вниз и по "результирующей" диагонали (см. рисунок выше).

## Корректировка диаграммы (для опытных пользователей).
```
Не закрывайте и не перезапускайте интерактивную консоль, пока не просмотрели полученную диаграмму!
```

Если после появления `>>>` консоль IDLE всё ещё открыта и не перезапущена вами (`CTRL+F6`), созданные в процессе работы программы объекты хранятся в оперативной памяти. Этим можно воспользоваться, пересоздав диаграмму без повторного выполнения целого скрипта. Если же запускать скрипт по новой, то придётся ждать, когда все значения LD заново сосчитаются, что может длиться очень долго. Инструкция полезна только в том случае, если на выходе - одна диаграмма. Сразу несколько так подправить не удастся.

Если вы знаете основы Python, то без труда разберётесь в объекте, формирующим диаграмму. Выведите этот объект на экран.
```
ld_heatmap
```

(дважды кликните на `squeezed text`).

К элементам этой структуры можно обращаться как к элементам обычных словарей и списков. Можно добавлять новые элементы и заменять имеющиеся. Разве что `del` не работает, но эту проблему можно обойти. Если вы - не Python-программист и не разбираетесь в структурах данных, то просто копируйте в интерактивную консоль приведённые ниже строки кода.

Всплывающие тексты в статье и презентации не отобразить — в этом случае максимальную информативность диаграммы можно обеспечить лишь надписями внутри квадратиков и вдоль осей. Тепловые карты с текстом, размещённым, как минимум, внутри квадратиков, далее будем называть аннотированными. Самая частая неприятность в таких диаграммах - когда надписи наползают друг на друга. Если диаграмма не очень большая (в пределах 50х50), то проблема решается уменьшением шрифта. В Plotly размер шрифта по умолчанию - 12, поэтому в качестве примера уменьшенного шрифта здесь и далее привожу 10. Смело пробуйте и другие значения. Если что, объект диаграммы можно будет пересоздать (см. ниже), тем самым сбросив ваши изменения.

В первую очередь, по мере увеличения тепловой карты начинают пересекаться LD-значения внутри квадратиков. Сделаем их поменьше.
```
for ann_num in range(len(ld_heatmap['layout']['annotations'])):
	ld_heatmap['layout']['annotations'][ann_num]['font']['size'] = 10
```

Кстати, есть альтернативный способ обращения к подструктурам Plotly-объектов - точечная запись.
```
for ann_num in range(len(ld_heatmap.layout.annotations)):
	ld_heatmap.layout.annotations[ann_num].font.size = 10
```

После каждого внесённого изменения предпросматривайте тепловую карту (при этом будет сохраняться только временный файл).
```
py.offline.plot(ld_heatmap)
```

Если размер матрицы ещё больше, то могут сливаться лейблы осей. Компактизируем их.
```
ld_heatmap['layout']['xaxis']['tickfont'] = {'size': 10}
```
```
ld_heatmap['layout']['yaxis']['tickfont'] = {'size': 10}
```

Кое-как ещё умещаются лейблы осей, но совсем не влезают надписи в квадратиках? Тогда последние удалим. Как уже говорилось ранее, с помощью `del` это не произвести, но можно схитрить - вместо ненужного объекта создать пустой объект того же типа. Взяв на вооружение этот способ, упраздняем список с настройками текстового заполнения квадратиков, тем самым освобождая квадратики от надписей.
```
ld_heatmap['layout']['annotations'] = []
```

Если матрица настолько крупная, что лейблы осей тоже стали выходить из берегов, то удаляем и их.
```
ld_heatmap['layout']['xaxis']['showticklabels'] = False
```
```
ld_heatmap['layout']['yaxis']['showticklabels'] = False
```

Получится тепловая карта такого вида:
![Пример крупной тепловой карты без надписей](https://raw.githubusercontent.com/PlatonB/ld-tools/master/LD_triangle_big_heatmap.png)

Опять что-то кажется кривым? Продолжайте описанные выше эксперименты с параметрами. Если же проделанная работа устраивает, можете сохранить обновлённую диаграмму в ту папку, путь к которой указывали в рамках пользовательского диалога. Учтите, что старая диаграмма при этом заменится на новую.
```
py.offline.plot(ld_heatmap, filename = html_file_path, auto_open=False)
```

Может так сложиться, что ваши эксперименты зашли слишком далеко и нужно вернуть тепловой карте тот вид, который она имела непосредственно после завершения работы программы. Ну или вы хотите, стартуя с нуля, вновь кастомизировать её до неузнаваемости. Тогда создадим объект диаграммы заново.

Объект аннотированной тепловой карты. Специальная функция из "фабрики фигур" Plotly сделает всю сложную низкоуровневую работу по построению этого объекта за нас.
```
ld_heatmap = ff.create_annotated_heatmap(ld_two_dim, x=rs_ids_srtd, y=rs_ids_srtd, hovertext = info_two_dim, hoverinfo = 'text', xgap=1, ygap=1, colorscale=color_map, reversescale=True)
```

Объект обычной тепловой карты. Соберём его вручную. Должен получиться словарь, состоящий из двух словареподобных объектов: один — с данными и настройками их отображения (`trace`), другой — с настройками не связанных с данными элементов диаграммы (`layout`). В нашем случае trace будет содержать значения LD, соответствующие всплывающие аннотации и некоторые кастомизации дизайна, а layout — параметры осей.
```
trace = go.Heatmap(z=ld_two_dim, hovertext=info_two_dim, hoverinfo='text', xgap=1, ygap=1, colorscale=color_map, reversescale=True, showscale=False)
```
```
layout = {'xaxis': {'showticklabels': False}, 'yaxis': {'showticklabels': False}}
```
```
ld_heatmap = {'data': [trace], 'layout': layout}
```

Перевернём тепловую карту по оси Y, чтобы она выглядела, как результат работы LDmatrix.
```
ld_heatmap['layout']['yaxis']['autorange'] = 'reversed'
```

Мы вернули диаграмме дефолтный вид. Посмотрим, что получилось.
```
py.offline.plot(ld_heatmap)
```

Обычная тепловая карта:
![Пример обычной тепловой карты до кастомизации](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_triangle_simple_heatmap_before_customization.png)

Аннотированная тепловая карта:
![Пример аннотированной тепловой карты до кастомизации](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_triangle_annotated_heatmap_before_customization.png)

Пересоздавая объект тепловой карты, поиграйте заодно с аргументами соответствующей функции. К примеру, настройте ширину промежутков между квадратиками (`xgap=2`, `ygap=2`), поменяйте цветовую схему (`colorscale='Hot'`).

Для поддержания схожести с диаграммами LDmatrix можно сделать наш вариант полностью квадратным. Для этого установим соотношение расстояний между делениями по X и делениями по Y 1:1, а также приведём к квадратной форме область построения (допустим, 1000х1000 пикселей).
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

Тогда созданная по новой диаграмма будет выглядеть так:
![Пример тепловой карты 1000х1000 пикселей с квадратными элементами](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_triangle_heatmap_after_customization.png)

Я охватил только самые напрашивающиеся идеи редактирования получаемых программой тепловых карт. Если вы обладаете колоссальным терпением и непреодолимым желанием производить тонкую настройку диаграмм Plotly, то эта ссылка для вас:
https://plot.ly/python/reference/
