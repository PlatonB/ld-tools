# Синопсис.
## Преимущества.
Основанное на собственном опыте сравнение _LD_tools_ с такими известными решениями, как _LDlink_, _Ensembl REST API_ и 
_Haploreg_.
- Работа с исходными файлами:
	- не нужно вручную загружать по одному файлу, программа сама сканирует всю папку с данными;
	- автоматически распознаётся столбец с идентификаторами SNP каждой таблицы.
- После подготовки кэша работу можно проводить полностью оффлайн. Никаких тормознутых API и 503-сюрпризов.
- Возможность тонкой настройки с помощью дружественного интерактивного диалога. Не требуется командная строка.
- Геномные координаты - только GRCh38. Нет характерной для биоинформатических тулзов путаницы по этому вопросу.
- Код на чистом Python с человеко-понятными именами переменных и подробными комментариями. Делаю всё, чтобы проект пригождался в том числе и для изучения программирования в биоинформатике.

## Перед началом работы.
1. Установите [IDLE](https://github.com/PlatonB/bioinformatic-python-scripts#установка-среды-разработки) или другую удобную вам среду разработки, поддерживающую Python 3. Либо можете, не используя среду разработки, запускать Python-программы из [терминала](https://github.com/PlatonB/ngs-pipelines#преодолеваем-страх-командной-строки-linux).
2. Разрешите [зависимости](https://github.com/PlatonB/ld-tools#установка-сторонних-компонентов).
3. Скачайте архив с программой:
```
Clone or download
```
(зелёная кнопка наверху страницы репозитория).

4. Распакуйте его в любую папку.

### Установка сторонних компонентов.
#### MongoDB.
Советую вначале ознакомиться с [основами работы в линуксовым терминале](https://github.com/PlatonB/ngs-pipelines#преодолеваем-страх-командной-строки-linux). Впрочем, если совсем лень, можете просто копировать, вставлять и запускать приведённые ниже команды. После установки настоятельно рекомендую перезагрузиться.

##### Ubuntu Linux.
([elementary OS](https://elementary.io/ru/)/KDE neon/Linux Mint)

Подключение официального репозитория MongoDB.
```
wget -qO - https://www.mongodb.org/static/pgp/server-4.2.asc | sudo apt-key add -
```
```
echo "deb [ arch=amd64 ] https://repo.mongodb.org/apt/ubuntu bionic/mongodb-org/4.2 multiverse" | sudo tee /etc/apt/sources.list.d/mongodb-org-4.2.list
```

Обновление индекса пакетов ОС.
```
sudo apt update
```

Собственно, установка MongoDB.
```
sudo apt install -y mongodb-org
```

Перманентный запуск MongoDB. Лучше так сделать, если планируете использовать _ld-tools_ и [high-perf-bio](https://github.com/PlatonB/high-perf-bio) часто.
```
systemctl enable mongod.service
```

Если вам не нужно эксплуатировать MongoDB-решения каждый день, то рекомендую команду, активирующую MongoDB до момента перезагрузки.
```
sudo service mongod start
```

##### Fedora Linux.
TBD.

#### Python-библиотеки.
Установка с помощью _pip_:
```
pip3 install pymongo plotly numpy --user
```

Установка с помощью [Conda](https://github.com/PlatonB/ngs-pipelines#установка-conda):
```
conda install pymongo plotly numpy
```

#### Примечание по поводу Windows.
Теоретически, после установки MongoDB и whl-пакетов _pymongo_, _plotly_ и _numpy_, программа должна работать. Но у меня сейчас Windows нет, и я пока не проверял. Надеюсь, кто-нибудь поделится опытом в [Issues](https://github.com/PlatonB/ld-tools/issues).

## Запуск.
1. Если пользуетесь _IDLE_, откройте в нём нужный компонент программы и запустите его кнопкой `F5`.
2. Введите [опции](https://github.com/PlatonB/bioinformatic-python-scripts#диалог-с-пользователем) в рамках интерактивного диалога.

В _IDLE_, когда работа программы завершится, появятся символы `>>>`.

## Важно!
- Мультиаллельные SNP не поддерживаются, поэтому будут проигнорированы любым компонентом этой программы.
- Подготовка данных, необходимых для быстрого оффлайн-вычисления LD, может занимать длительное время и сопровождаться небольшими лагами ОС. Но выполняться она будет только при первом запуске программы. Не отключайте интернет, и просто перетерпите сутки-двое:).
- По любым вопросам пишите в [Issues](https://github.com/PlatonB/ld-tools/issues) (потребуется регистрация на Github).
- Другие полезные биоинформатические инструменты моей разработки можете найти в репозитории https://github.com/PlatonB/bioinformatic-python-scripts.

# LD_area.
Выводит для каждого из запрашиваемых SNP все SNP в пределах двух порогов: не ниже определённого LD и не дальше заданных фланков.

## Ключевые особенности.
- Файлы, получающиеся на выходе (по выбору пользователя):
	- JSON. Список словарей со значениями LD и физической дистанции, а также идентификаторами и наиболее важными характеристики каждого SNP. Прекрасный вариант, если хотите обрабатывать результаты скриптами.
![Пример JSON](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_area_JSON_output.png)
	- Набор одних лишь refSNPID. Такая разновидность вывода лучше всего подходит для дальнейшего аннотирования специализированными программами (_Ensembl VEP_ и др.).
	- TSV. Для ознакомления с результатами, редактирования в офисных пакетах и скриптовой обработки.
- Распределение результатов по хромосомным подпапкам.

# LD_triangle.
Строит матрицы значений неравновесия по сцеплению для пар SNP. Каждая матрица создаётся по всем парным сочетаниям соответствующего единого набора refSNPID. Важно отметить, что если уже найдено LD для пары SNPx-SNPy, то обработка пары SNPy-SNPx допущена не будет. Поэтому вместо квадратов с результатами получаются прямоугольные треугольники. Программа является написанной с нуля десктопной реализацией идеи компонента _LDmatrix_ веб-сервиса _LDlink_.

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

- (ИМХО!) Более наглядные и симпатичные диаграммы (можете оспорить это в [Issues](https://github.com/PlatonB/ld-tools/issues):smile:).

## Подходящие цветовые палитры.
Подобрал такие палитры _Plotly_, с которыми LD-карты не будут смотреться вырвиглазно.

_algae, amp, blues, blugrn, bluyl, brwnyl, bugn, bupu, burg, burgyl, darkmint, deep, dense, emrld, gnbu, greens, greys, magenta, matter, mint, oranges, orrd, oryel, peach, pinkyl, pubu, pubugn, purd, purp, purples, purpor, rdpu, redor, reds, speed, sunset, sunsetdark, teal, tealgrn, tempo, turbid, ylgn, ylgnbu, ylorbr, ylorrd_

## Результаты.
На пересечении строки с одним SNP и столбца с другим отображается значение LD. Оно может быть представлено в виде вещественного числа, степени насыщенности окраски квадратика, либо сочетанием одного и другого.
![Пример тепловой карты с указанием закономерностей физического расстояния](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_triangle_annotated_heatmap_with_explanations.png)
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

Всплывающие тексты в статье и презентации не отобразить — в этом случае максимальную информативность диаграммы можно обеспечить лишь надписями внутри квадратиков и вдоль осей. Тепловые карты с текстом, размещённым, как минимум, внутри квадратиков, далее будем называть аннотированными. Самая частая неприятность в таких диаграммах - когда надписи наползают друг на друга. Если диаграмма не очень большая (в пределах 50х50), то проблема решается уменьшением шрифта. В _Plotly_ размер шрифта по умолчанию - 12, поэтому в качестве примера уменьшенного шрифта здесь и далее привожу 10. Смело пробуйте и другие значения. Если что, объект диаграммы можно будет пересоздать (см. ниже), тем самым сбросив ваши изменения.

В первую очередь, по мере увеличения тепловой карты начинают пересекаться LD-значения внутри квадратиков. Сделаем их поменьше.
```
for ann_num in range(len(ld_heatmap['layout']['annotations'])):
	ld_heatmap['layout']['annotations'][ann_num]['font']['size'] = 10
```

Кстати, есть альтернативный способ обращения к подструктурам _Plotly_-объектов - точечная запись.
```
for ann_num in range(len(ld_heatmap.layout.annotations)):
	ld_heatmap.layout.annotations[ann_num].font.size = 10
```

После каждого внесённого изменения предпросматривайте тепловую карту.
```
ld_heatmap.show()
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
![Пример крупной тепловой карты без надписей](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_triangle_big_heatmap.png)

Опять что-то кажется кривым? Продолжайте описанные выше эксперименты с параметрами. Если же проделанная работа устраивает, можете сохранить обновлённую диаграмму в ту папку, путь к которой указывали в рамках пользовательского диалога. Учтите, что старая диаграмма при этом заменится на новую.
```
ld_heatmap.write_html(html_file_path)
```

Может так сложиться, что ваши эксперименты зашли слишком далеко и нужно вернуть тепловой карте тот вид, который она имела непосредственно после завершения работы программы. Ну или вы хотите, стартуя с нуля, вновь кастомизировать её до неузнаваемости. Тогда создадим объект диаграммы заново.

Объект аннотированной тепловой карты. Специальная функция из "фабрики фигур" _Plotly_ сделает всю сложную низкоуровневую работу по построению этого объекта за нас.
```
ld_heatmap = ff.create_annotated_heatmap(ld_two_dim, x=rs_ids_srtd, y=rs_ids_srtd, hovertext=info_two_dim, hoverinfo='text', xgap=1, ygap=1, colorscale=color_pal)
```

Объект обычной тепловой карты. Соберём его вручную. Должен получиться словареподобный объект, состоящий из двух объектов: один — с данными и настройками их отображения (`trace`), другой — с настройками не связанных с данными элементов диаграммы (`layout`). В нашем случае trace будет содержать значения LD, соответствующие всплывающие аннотации и некоторые кастомизации дизайна, а layout — параметры осей.
```
trace = go.Heatmap(z=ld_two_dim, hovertext=info_two_dim, hoverinfo='text', xgap=1, ygap=1, colorscale=color_pal, showscale=False)
```
```
layout = go.Layout(xaxis={'showticklabels': False}, yaxis={'showticklabels': False})
```
```
ld_heatmap = go.Figure(data=[trace], layout=layout)
```

Перевернём тепловую карту по оси Y, чтобы она выглядела, как результат работы _LDmatrix_.
```
ld_heatmap['layout']['yaxis']['autorange'] = 'reversed'
```

Мы вернули диаграмме дефолтный вид. Посмотрим, что получилось.
```
ld_heatmap.show()
```

Обычная тепловая карта:
![Пример обычной тепловой карты до кастомизации](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_triangle_simple_heatmap_before_customization.png)

Аннотированная тепловая карта:
![Пример аннотированной тепловой карты до кастомизации](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_triangle_annotated_heatmap_before_customization.png)

Пересоздавая объект тепловой карты, поиграйте заодно с аргументами соответствующей функции. К примеру, уберите промежутки между квадратиками (`xgap=0`, `ygap=0`), поменяйте цветовую палитру (`colorscale='amp'`).

Для поддержания схожести с диаграммами _LDmatrix_ можно сделать наш вариант полностью квадратным. Для этого установим соотношение расстояний между делениями по X и делениями по Y 1:1, а также приведём к квадратной форме область построения (допустим, 1000х1000 пикселей).
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

Объект квадратной тепловой карты, в целом, готов. Но четвёртая версия _Plotly_ порождает визуальную шероховатость: серый бэкграунд будет чуточку вылезать за бока диаграммы. Выход из положения — сделать бэкграунд белым.
```
ld_heatmap['layout']['plot_bgcolor'] = 'rgba(0,0,0,0)'
```

Созданная по новой диаграмма будет выглядеть так:
![Пример тепловой карты 1000х1000 пикселей с квадратными элементами](https://raw.githubusercontent.com/PlatonB/ld-tools/master/gallery/LD_triangle_heatmap_after_customization.png)

Я охватил только самые напрашивающиеся идеи редактирования получаемых программой тепловых карт. Если вы обладаете колоссальным терпением и непреодолимым желанием производить тонкую настройку диаграмм _Plotly_, то эта ссылка для вас:
https://plot.ly/python/reference/
