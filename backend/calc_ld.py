__version__ = 'V5.0'

def calc_ld(snp_1_genotypes, snp_2_genotypes):
        '''
        Вычисление значений
        LD (r2 и D') для двух
        биаллельных SNP по генотипам.
        Функция запрашивает на вход
        лишь два списка или кортежа
        генотипов, что делает калькулятор
        легко встраиваемым в другие
        биоинформатические программы.
        Предполагается, что генотипы
        были ранее отобраны по
        принадлежности соответствующих
        сэмплов полу и популяциям.
        Генотипы в этих массивах
        должны быть исключительно
        исключительно в одиночном виде:
        вместо 1|1, 1|0 - 1, 1, 1, 0.
        0 означает референсный
        аллель, а 1 - альтернативный.
        '''
        
        #Определение общего количества гаплотипов.
        #Оно равно общему количеству отдельных генотипов.
        #Это значение константно в пределах
        #каждой таблицы 1000 Gemomes.
        htypes = list(zip(snp_1_genotypes, snp_2_genotypes))
        
        #Рассчёт количества гаплотипов
        #alt+alt и каждого аллеля.
        htype_alt_alt_quan = htypes.count((1, 1))
        snp_1_alt_quan = snp_1_genotypes.count(1)
        snp_1_ref_quan = snp_1_genotypes.count(0)
        snp_2_alt_quan = snp_2_genotypes.count(1)
        snp_2_ref_quan = snp_2_genotypes.count(0)
        
        #Рассчёт количества и
        #частоты гаплотипа alt+alt
        #и частоты каждого аллеля.
        htype_quan = len(htypes)
        alt_alt_htype_freq = htype_alt_alt_quan / htype_quan
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
        #SNP, т.е. у которого все генотипы
        #одинаковы: только 1|1 или только 0|0.
        #Расчёт D' в этом случае невозможен
        #из-за попадания нуля в знаменатель.
        #Поскольку None-значения пока вызывают
        #проблемы на уровне фронтендов, при обработке
        #мономорфных SNP будет возвращён D', равный 0.
        if d >= 0:
                d_max = min(snp_1_alt_freq * snp_2_ref_freq,
                            snp_1_ref_freq * snp_2_alt_freq)
                try:
                        d_prime = d / d_max
                except ZeroDivisionError:
                        d_prime = 0
        else:
                d_min = max(-snp_1_alt_freq * snp_2_alt_freq,
                            -snp_1_ref_freq * snp_2_ref_freq)
                try:
                        d_prime = d / d_min
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
                r_square = (d ** 2) / (snp_1_alt_freq * snp_1_ref_freq *
                                       snp_2_alt_freq * snp_2_ref_freq)
        else:
                r_square = 0
                
        #Оба значения LD и частоты
        #альтернативных аллелей
        #организуем в словарь.
        results = {'r_square': round(r_square, 4),
                   'd_prime': round(d_prime, 4),
                   'snp_1_alt_freq': round(snp_1_alt_freq, 4),
                   'snp_2_alt_freq': round(snp_2_alt_freq, 4)}
        
        return results
