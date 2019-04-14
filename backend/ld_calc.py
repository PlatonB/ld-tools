__version__ = 'V3.0'

import os, json

def ld_calc(snp_1_row, snp_2_row, sample_indices):
        '''
        Вычисление значений LD (r2 и D') для
        двух SNP по фазированным генотипам.
        '''
        
        #Получение для выбранных сэмплов фазированных генотипов,
        #соответствующих каждому из двух биаллельных SNP.
        #По ним можно понять, альтернативное или
        #референсное основание находится в определённом
        #локусе у индивида из выборки 1000 Genomes.
        #Фазированные генотипы бывают такими: 1|1, 1|0, 0|1, 0|0, 1, 0.
        #То, что в каждом, если не рассматривать X и Y-хромосомы,
        #два значения - следствие диплоидности у человека.
        #К примеру, 0|1 показывает, что в данной позиции у
        #этого индивида на одной копии хромосомы - референсный
        #аллель (0), а на другой - альтернативный (1).
        snp_1_phased_genotypes = [snp_1_row[sample_index] for sample_index in sample_indices]
        snp_2_phased_genotypes = [snp_2_row[sample_index] for sample_index in sample_indices]
        
        #Разделение фазированных генотипов на генотипы каждой копии хромосомы.
        #Для хромосом X (♂) и Y генотипы были и останутся одиночными.
        snp_1_genotypes, snp_2_genotypes = [], []
        for phased_genotype_num in range(len(snp_1_phased_genotypes)):
                snp_1_genotypes += snp_1_phased_genotypes[phased_genotype_num].split('|')
                snp_2_genotypes += snp_2_phased_genotypes[phased_genotype_num].split('|')
                
        #Определение общего количества гаплотипов.
        #Оно равно общему количеству отдельных генотипов.
        #Это значение константно в пределах
        #каждой таблицы 1000 Gemomes.
        htype_quan = len(snp_1_genotypes)
        
        #Рассчёт количества каждого гаплотипа и каждого аллеля.
        htype_alt_alt_quan, htype_alt_ref_quan, htype_ref_alt_quan, htype_ref_ref_quan = 0, 0, 0, 0
        snp_1_alt_quan, snp_1_ref_quan, snp_2_alt_quan, snp_2_ref_quan = 0, 0, 0, 0
        for htype in zip(snp_1_genotypes, snp_2_genotypes):
                if htype[0] == '1' and htype[1] == '1':
                        htype_alt_alt_quan += 1
                        snp_1_alt_quan += 1
                        snp_2_alt_quan += 1
                elif htype[0] == '1' and htype[1] == '0':
                        htype_alt_ref_quan += 1
                        snp_1_alt_quan += 1
                        snp_2_ref_quan += 1
                elif htype[0] == '0' and htype[1] == '1':
                        htype_ref_alt_quan += 1
                        snp_1_ref_quan += 1
                        snp_2_alt_quan += 1
                elif htype[0] == '0' and htype[1] == '0':
                        htype_ref_ref_quan += 1
                        snp_1_ref_quan += 1
                        snp_2_ref_quan += 1
                        
        #Рассчёт частоты каждого гаплотипа и каждого аллеля.
        alt_alt_htype_freq = htype_alt_alt_quan / htype_quan
        alt_ref_htype_freq = htype_alt_ref_quan / htype_quan
        ref_alt_htype_freq = htype_ref_alt_quan / htype_quan
        ref_ref_htype_freq = htype_ref_ref_quan / htype_quan
        snp_1_alt_freq = snp_1_alt_quan / htype_quan
        snp_1_ref_freq = snp_1_ref_quan / htype_quan
        snp_2_alt_freq = snp_2_alt_quan / htype_quan
        snp_2_ref_freq = snp_2_ref_quan / htype_quan
        
        #D - стандартное значение LD.
        #Если оно равно нулю, то это
        #означает равновесие по сцеплению.
        #Неравенство нулю будет говорить
        #о неравновесии по сцеплению.
        #D может быть от -0.25 до 0.25.
        d = alt_alt_htype_freq - snp_1_alt_freq * snp_2_alt_freq
        
        #D' - нормализованное значение LD.
        #Частота одного из аллелей в редких
        #случаях может быть равной нулю.
        #Это означает, что попался мономорфный
        #SNP, т.е. у которого все генопипы
        #одинаковы: только 1|1 или только 0|0.
        #Расчёт D' в этом случае невозможен
        #из-за попадания нуля в знаменатель.
        #Поскольку None-значения пока вызывают
        #проблемы на уровне фронтендов, при обработке
        #мономорфных SNP будет возвращён D', равный 0.
        if d >= 0:
                d_max = min(snp_1_alt_freq * snp_2_ref_freq, snp_1_ref_freq * snp_2_alt_freq)
                try:
                        d_prime = round(d / d_max, 5)
                except ZeroDivisionError:
                        d_prime = 0
        else:
                d_min = max(-snp_1_alt_freq * snp_2_alt_freq, -snp_1_ref_freq * snp_2_ref_freq)
                try:
                        d_prime = round(d / d_min, 5)
                except ZeroDivisionError:
                        d_prime = 0
                        
        #r2 - значение LD в виде квадрата
        #коэффициента корреляции Пирсона.
        #Если r2 = 0, то это полное равновесие по сцеплению,
        #если r2 = 1, то тогда полное неравновесие.
        #Если один из SNP - мономорфный, то
        #в качестве значения r2 возвратится 0.
        #Эта ситуация подробно описана выше для D'.
        if d_prime != 0:
                r_square = round((d ** 2) / (snp_1_alt_freq * snp_1_ref_freq * snp_2_alt_freq * snp_2_ref_freq), 5)
        else:
                r_square = 0

        #Оба значения LD, а также
        #частоты альтернативных аллелей,
        #организуем в словарь и возвращаем.
        results = {'r_square': r_square,
                   'd_prime': d_prime,
                   'snp_1_alt_freq': round(snp_1_alt_freq, 5),
                   'snp_2_alt_freq': round(snp_2_alt_freq, 5)}
        return results
