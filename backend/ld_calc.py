__version__ = 'V1.0'

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
        #К примеру, 0|1 показывает, что в данной позиции
        #у этого индивида на одной копии хромосомы -
        #референсный аллель, а на другой - альтернативный.
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
        haplotype_quan = len(snp_1_genotypes)
        
        #Рассчёт количества каждого гаплотипа и каждого аллеля.
        haplotype_11_quan, haplotype_12_quan, haplotype_21_quan, haplotype_22_quan = 0, 0, 0, 0
        a1_allele_quan, a2_allele_quan, b1_allele_quan, b2_allele_quan = 0, 0, 0, 0
        for haplotype in zip(snp_1_genotypes, snp_2_genotypes):
                if haplotype[0] == '1' and haplotype[1] == '1':
                        haplotype_11_quan += 1
                        a1_allele_quan += 1
                        b1_allele_quan += 1
                elif haplotype[0] == '1' and haplotype[1] == '0':
                        haplotype_12_quan += 1
                        a1_allele_quan += 1
                        b2_allele_quan += 1
                elif haplotype[0] == '0' and haplotype[1] == '1':
                        haplotype_21_quan += 1
                        a2_allele_quan += 1
                        b1_allele_quan += 1
                elif haplotype[0] == '0' and haplotype[1] == '0':
                        haplotype_22_quan += 1
                        a2_allele_quan += 1
                        b2_allele_quan += 1
                        
        #Рассчёт частоты каждого гаплотипа и каждого аллеля.
        x11 = haplotype_11_quan / haplotype_quan
        x12 = haplotype_12_quan / haplotype_quan
        x21 = haplotype_21_quan / haplotype_quan
        x22 = haplotype_22_quan / haplotype_quan
        p1 = a1_allele_quan / haplotype_quan
        p2 = a2_allele_quan / haplotype_quan
        q1 = b1_allele_quan / haplotype_quan
        q2 = b2_allele_quan / haplotype_quan
        
        #D - стандартное значение LD.
        #Если оно равно нулю, то это
        #означает равновесие по сцеплению.
        #Неравенство нулю будет говорить
        #о неравновесии по сцеплению.
        #D может быть от -0.25 до 0.25.
        d = x11 - p1 * q1
        
        #D' - нормализованное значение LD.
        if d > 0:
                d_max = min(p1 * q2, p2 * q1)
                d_prime = d / d_max
        elif d < 0:
                d_min = max(-p1 * q1, -p2 * q2)
                d_prime = d / d_min
        else:
                d_prime = 0
                
        #r2 - значение LD в виде квадрата
        #коэффициента корреляции Пирсона.
        #Если r2 = 0, то это полное равновесие по сцеплению,
        #если r2 = 1, то тогда полное неравновесие.
        if d_prime != 0:
                r_square = (d ** 2) / (p1 * p2 * q1 * q2)
        else:
                r_square = 0

        #Оба значения LD организуем в словарь и возвращаем.
        ld_vals = {'r_square': r_square, 'd_prime': d_prime}
        return ld_vals
