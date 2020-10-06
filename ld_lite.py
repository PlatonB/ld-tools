__version__ = 'V1.0'

class NotRsIdError(Exception):
        '''
        Несоответствие идентификатора
        варианта стандарту reference SNP ID.
        '''
        def __init__(self, rs_id):
                err_msg = f'{rs_id} is non-rs identifier'
                super().__init__(err_msg)
                
class NotInIntgenConvDbError(Exception):
        '''
        Невхождение варианта в основанную
        на данных 1000 Genomes базу
        конвертации rsIDs и координат.
        '''
        def __init__(self, rs_id):
                err_msg = f'{rs_id} is not available in 1000 Genomes'
                super().__init__(err_msg)
                
class DifChrsError(Exception):
        '''
        Если изучаемая пара вариантов
        находится на разных хромосомах,
        для неё невозможно будет
        сосчитать LD и дистанцию.
        '''
        def __init__(self, rs_id_1, rs_id_2):
                err_msg = f'{rs_id_1} and {rs_id_2} belong to different chromosomes'
                super().__init__(err_msg)
                
def check_rs_id(rs_id, cursor):
        '''
        Проверка исходного варианта
        на валидность, получение
        его хромосомы и позиции.
        '''
        if re.search(r'rs\d+', rs_id) is None:
                raise NotRsIdError(rs_id)
        cursor.execute(f'SELECT CHROM, POS FROM variants WHERE ID = "{rs_id}"')
        var_basic_info = cursor.fetchone()
        if var_basic_info is None:
                raise NotInIntgenConvDbError(rs_id)
        return var_basic_info

####################################################################################################

import sys, locale, os, re, sqlite3
from tabulate import tabulate

#Подавление формирования питоновского кэша с
#целью предотвращения искажения результатов.
sys.dont_write_bytecode = True

from cli.ld_lite_cli_ru import add_args_ru
from cli.ld_lite_cli_en import add_args_en
from backend.prep_intgen_data import prep_intgen_data
from backend.get_sample_names import get_sample_names
from pysam import VariantFile
from backend.calc_ld import calc_ld

#Обработка аргументов командной строки.
if locale.getdefaultlocale()[0][:2] == 'ru':
        args = add_args_ru(__version__)
else:
        args = add_args_en(__version__)
        
intgen_dir_path = os.path.normpath(args.intgen_dir_path)
if args.skip_intgen_data_ver:
        intgen_convdb_path = os.path.join(intgen_dir_path,
                                          'conversion.db')
else:
        intgen_convdb_path = prep_intgen_data(intgen_dir_path)
if args.gend_names == 'male':
        gend_names = ('male',)
elif args.gend_names == 'female':
        gend_names = ('female',)
else:
        gend_names = ('male', 'female')
pop_names = tuple(args.pop_names.upper().split(','))
sample_names = get_sample_names(gend_names,
                                pop_names,
                                intgen_convdb_path)

#Подключаемся к базе, содержащей имена хромосом,
#позиции и rsIDs биаллельных вариантов из 1000
#Genomes. Проверяем на корректность два исходных
#варианта, определяем их локализацию и убеждаемся
#в том, что они оба принадлежат одной хромосоме.
with sqlite3.connect(intgen_convdb_path) as conn:
        cursor = conn.cursor()
        var_1_basic_info = check_rs_id(args.rs_id_1, cursor)
        var_2_basic_info = check_rs_id(args.rs_id_2, cursor)
        cursor.close()
if var_1_basic_info[0] != var_2_basic_info[0]:
        raise DifChrsError(args.rs_id_1, args.rs_id_2)
chrom, var_1_pos = var_1_basic_info
var_2_pos = var_2_basic_info[1]

#Поскольку известно, на какой хромосоме
#лежит поступившая на вход пара вариантов,
#можно сразу открыть архив именно с
#относящимися к этой хромосоме данными
#1000 Genomes. Из него извлечём генотипы
#и некоторые аннотации исходных вариантов.
#Тот ли вариант нашёлся в 1000 Genomes,
#дополнительно проверяется по rsID.
with VariantFile(os.path.join(intgen_dir_path,
                              f'{chrom}.vcf.gz')) as intgen_vcf_opened:
        var_1_genotypes, var_2_genotypes = [], []
        for intgen_rec in intgen_vcf_opened.fetch(chrom,
                                                  var_1_pos - 1,
                                                  var_1_pos):
                if intgen_rec.id != args.rs_id_1:
                        continue
                var_1_alleles = intgen_rec.ref + '/' + intgen_rec.alts[0]
                for sample_name in sample_names:
                        try:
                                var_1_genotypes += intgen_rec.samples[sample_name]['GT']
                        except KeyError:
                                continue
                break
        for intgen_rec in intgen_vcf_opened.fetch(chrom,
                                                  var_2_pos - 1,
                                                  var_2_pos):
                if intgen_rec.id != args.rs_id_2:
                        continue
                var_2_alleles = intgen_rec.ref + '/' + intgen_rec.alts[0]
                for sample_name in sample_names:
                        try:
                                var_2_genotypes += intgen_rec.samples[sample_name]['GT']
                        except KeyError:
                                continue
                break
        
#Получение значений LD с помощью
#оффлайн-калькулятора. Полезный
#побочный продукт функции - частоты
#альтернативного аллеля вариантов.
trg_vals = calc_ld(var_1_genotypes,
                   var_2_genotypes)

#Вывод результатов, прилично оформленных
#с помощью конструктора таблиц tabulate.
print(tabulate([['chrom', chrom, chrom],
                ['hg38_pos', var_1_pos, var_2_pos],
                ['alleles', var_1_alleles, var_2_alleles],
                ['alt_freq', trg_vals['var_1_alt_freq'], trg_vals['var_2_alt_freq']]],
               headers=[tabulate([['r2', trg_vals['r_square']],
                                  ["D'", trg_vals['d_prime']],
                                  ['abs_dist', abs(var_1_pos - var_2_pos)]],
                                 tablefmt='fancy_grid',
                                 disable_numparse=True),
                        f'\n\n\n{args.rs_id_1}', f'\n\n\n{args.rs_id_2}'],
               tablefmt='fancy_grid'))
