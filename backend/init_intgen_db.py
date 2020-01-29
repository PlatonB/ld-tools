__version__ = 'V1.0'

from pymongo import MongoClient

#Создание Python-объектов
#клиента MongoDB, базы данных
#1000 Genomes и коллекции сэмплов.
client = MongoClient()
intgen_db_obj = client['1000_genomes_db']
intgen_sampcoll_obj = intgen_db_obj['samples']
