import numpy as np
import datetime
import sqlite3
import matplotlib.pyplot as plt
import csv
import prettytable
from sqlite3 import OperationalError
import sys, os
from collections import Counter


class DataProcessing:
    def __init__(self):
        # self.data_DF = np.loadtxt('./data/dengue_c_her2010.dat')
        # self.data_DHF = np.loadtxt('./data/dengue_h_her2010.dat')
        # self.dates_excel_DF = self.data_DF[:, 0]
        # self.dates_excel_DHF = self.data_DHF[:, 1]
        # self.day_zero_classic_excel = 40386
        # self.day_zero_hemorrhagic_excel = 40386
        self.result_DF = []
        self.result_DHF = []
        self.date_frecuency_DF = []
        self.date_frecuency_DHF = []
        self.frecuency_per_week_DF = []
        self.frecuency_per_week_DHF = []
        os.system("sh ./data/setup.sh")
        self.connection = sqlite3.connect('./data/dengue_data_2010.sqlite')
        self.cur = self.connection.cursor()

    @staticmethod
    def excel_date2python_datetime(excel_date):

        temp = datetime.datetime(1900, 1, 1)
        delta = datetime.timedelta(days=excel_date)
        return temp + delta

    def excel_data_dates2python_datetime(self, dates_excel):

        dates_py = []
        for day in dates_excel:
            date_py = self.excel_date2python_datetime(day)
            dates_py.append(date_py)
        return dates_py

    @staticmethod
    def indices(lst, element_to_match):
        result = []
        offset = -1
        while True:
            try:
                offset = lst.index(element_to_match, offset + 1)
            except ValueError:
                return result
            result.append(offset)

    def create_nominal_data(self):
        """
        :return: the nominal database in sql format.
        """
        try:
            self.cur.execute("CREATE TABLE nominal_data (\
                        OBJECTID, Join_Count, TARGET_FID, ID_HILLO, ID_SON,\
                        ANYO,\
                        FOLIO, SEXO, EDAD, D_H_, U_NOTIFICA, J_N_, DOMICILIO, \
                        COLONIA, LOCALIDAD, MUNICIPIO, I_CUADRO, I_CUADRO2, \
                        S_E__INICI, J_S_, HOSP, HEMORRAGIA, PLAQUETAS,\
                        DEFUNCION, FOLIO_DENG, DIAS_EVOLU,\
                        IgM_ELISA, IgG_ELISA, NS1, AISLAMIENT, FD_FHD,\
                        DIAGNOSTIC,\
                        OBSERVACIO, ID_HILLO_1, EDAD2, CUANTIFICA, GRAVE, \
                        CVEGEO,\
                        CVEGEO_1, POB1, POB26_R, POB27_R, POB30_R, ECO25_R,\
                        SALUD2_R, EDU34_R, MIG11_R, VIV2_R, VIV4_R, VIV9_R, \
                        VIV17_R, VIV26_R, VIV27_R, VIV28_R, VIV31_R, VIV32_R, \
                        VIV33_R, VIV34_R, VIV35_R, VIV36_R, GRADO_MARG, \
                        GRADO_MA_1, hectarea, DENSI_POB, VIV2, INDEX_NORM, \
                        INDEX_NO_1, x, y);")
        except OperationalError:
            None
        with open('./data/casos_dengue_2010_AGEB2010.csv', 'rb') as file_in:
            data_record_lines = csv.DictReader(file_in)
            to_db = [
                (
                    i["OBJECTID"], i["Join_Count"], i["TARGET_FID"],
                    i["ID_HILLO"], i["ID_SON"], i["ANYO"], i["FOLIO"],
                    i["SEXO"], i["EDAD"], i["D_H_"], i["U_NOTIFICA"],
                    i["J_N_"], i["DOMICILIO"], i["COLONIA"],
                    i["LOCALIDAD"],
                    i["MUNICIPIO"], i["I_CUADRO"], i["I_CUADRO2"],
                    i["S_E__INICI"], i["J_S_"], i["HOSP"],
                    i["HEMORRAGIA"],
                    i["PLAQUETAS"], i["DEFUNCION"], i["FOLIO_DENG"],
                    i["DIAS_EVOLU"], i["IgM_ELISA"], i["IgG_ELISA"],
                    i["NS1"],
                    i["AISLAMIENT"], i["FD_FHD"], i["DIAGNOSTIC"],
                    i["OBSERVACIO"], i["ID_HILLO_1"], i["EDAD2"],
                    i["CUANTIFICA"], i["GRAVE"], i["CVEGEO"],
                    i["CVEGEO_1"],
                    i["POB1"], i["POB26_R"], i["POB27_R"],
                    i["POB30_R"],
                    i["ECO25_R"], i["SALUD2_R"], i["EDU34_R"],
                    i["MIG11_R"],
                    i["VIV2_R"], i["VIV4_R"], i["VIV9_R"],
                    i["VIV17_R"],
                    i["VIV26_R"], i["VIV27_R"], i["VIV28_R"],
                    i["VIV31_R"],
                    i["VIV32_R"], i["VIV33_R"], i["VIV34_R"],
                    i["VIV35_R"],
                    i["VIV36_R"], i["GRADO_MARG"], i["GRADO_MA_1"],
                    i["hectarea"], i["DENSI_POB"], i["VIV2"],
                    i["INDEX_NORM"], i["INDEX_NO_1"], i["x"], i["y"]
                    )
                for i in data_record_lines
                ]
        self.cur.executemany("INSERT INTO nominal_data(\
                    OBJECTID, Join_Count, TARGET_FID, ID_HILLO, ID_SON,\
                    ANYO, FOLIO, SEXO, EDAD, D_H_, U_NOTIFICA, J_N_, \
                    DOMICILIO, COLONIA, LOCALIDAD, MUNICIPIO, I_CUADRO, \
                    I_CUADRO2, S_E__INICI, J_S_, HOSP, HEMORRAGIA, \
                    PLAQUETAS, DEFUNCION, FOLIO_DENG, DIAS_EVOLU,\
                    IgM_ELISA, IgG_ELISA, NS1, AISLAMIENT, FD_FHD, DIAGNOSTIC,\
                    OBSERVACIO, ID_HILLO_1, EDAD2, CUANTIFICA, GRAVE, CVEGEO,\
                    CVEGEO_1, POB1, POB26_R, POB27_R, POB30_R, ECO25_R, \
                    SALUD2_R, EDU34_R, MIG11_R, VIV2_R, VIV4_R, VIV9_R, \
                    VIV17_R, VIV26_R, VIV27_R, VIV28_R, VIV31_R, VIV32_R, \
                    VIV33_R, VIV34_R, VIV35_R, VIV36_R, GRADO_MARG,\
                    GRADO_MA_1, hectarea, DENSI_POB, VIV2, INDEX_NORM, \
                    INDEX_NO_1, x, y)\
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?,\
                            ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,\
                            ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,\
                            ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,\
                            ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,\
                            ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,\
                            ?, ?, ?, ?, ?, ?, ?, ?, ?);", to_db)
        self.connection.commit()
        return None

#    @staticmethod
    def create_table_incidence_data(self):
        """
        :return: the incidence case database in sql format.
        """
        try:
            self.cur.execute("CREATE TABLE incidence_data \
                        (case_date, case_week, fever, x, y);")
        except OperationalError:
            None

        with open('./data/casos_dengue_2010_AGEB2010.csv', 'rb') as file_in:
            data_record_lines = csv.DictReader(file_in)
            to_db = [
                (
                    i["I_CUADRO"], i["S_E__INICI"], i["FD_FHD"], i["x"], i["y"]
                    )
                for i in data_record_lines
                ]
        self.cur.executemany("INSERT INTO incidence_data(\
                    case_date, case_week, fever, x, y)\
                    VALUES (?, ?, ?, ?, ?);", to_db)
        self.connection.commit()
        return None

    def hemorrhagic_query_date_week(self):
        self.create_nominal_data()
        self.create_table_incidence_data()
        query_buffer_DF = "SELECT case_date, case_week, \
                            fever FROM incidence_data WHERE fever = 'FD'"
        self.cur.execute(query_buffer_DF)
        table_incidence = prettytable.from_db_cursor(self.cur)
        table_incidence.set_style(prettytable.PLAIN_COLUMNS)
        # print table_incidence
        table_incidence.header = True
        table_incidence.align = "l"
        table_incidence.left_padding_width = 0
        table_incidence.right_padding_width = 1
        string = table_incidence.get_string()
        string = str(string)
        result_DF = [tuple(filter(None, map(str.strip, split_line)))
                     for line in string.splitlines()
                     for split_line in [line.split(" ")]
                     if len(split_line) > 1]
        with open('./data/incidence_data_DF.csv', 'wb') as outcsv:
            writer = csv.writer(outcsv)
            writer.writerows(result_DF)
        #
        # frecuency per day counting
        #
        query_buffer_DHF = "SELECT case_date, case_week, \
                           fever FROM incidence_data WHERE fever = 'FHD'"
        self.cur.execute(query_buffer_DHF)
        table_incidence = prettytable.from_db_cursor(self.cur)
        table_incidence.set_style(prettytable.PLAIN_COLUMNS)
        # print table_incidence
        table_incidence.header = True
        table_incidence.align = "l"
        table_incidence.left_padding_width = 0
        table_incidence.right_padding_width = 1
        string = table_incidence.get_string()
        string = str(string)
        result_DHF = [tuple(filter(None, map(str.strip, splitline)))
                      for line in string.splitlines()
                      for splitline in [line.split(" ")]
                      if len(splitline) > 1]
        with open('./data/incidence_data_DHF.csv', 'wb') as outcsv:
            writer = csv.writer(outcsv)
            writer.writerows(result_DHF)

        self.result_DHF = result_DHF
        self.result_DF = result_DF
        return result_DF, result_DHF
    
    def hemorrhagic_query(self):
        self.create_nominal_data()
        self.create_table_incidence_data()
        query_buffer_DF = "SELECT case_date, \
                            fever FROM incidence_data WHERE fever = 'FD'"
        self.cur.execute(query_buffer_DF)
        table_incidence = prettytable.from_db_cursor(self.cur)
        table_incidence.set_style(prettytable.PLAIN_COLUMNS)
        # print table_incidence
        table_incidence.header = True
        table_incidence.align = "l"
        table_incidence.left_padding_width = 0
        table_incidence.right_padding_width = 1
        string = table_incidence.get_string()
        string = str(string)
        result_DF = [tuple(filter(None, map(str.strip, split_line)))
                     for line in string.splitlines()
                     for split_line in [line.split(" ")]
                     if len(split_line) > 1]
        with open('./data/incidence_data_DF.csv', 'wb') as outcsv:
            writer = csv.writer(outcsv)
            writer.writerows(result_DF)
        #
        # frecuency per day counting
        #
        query_buffer_DHF = "SELECT case_date, \
                           fever FROM incidence_data WHERE fever = 'FHD'"
        self.cur.execute(query_buffer_DHF)
        table_incidence = prettytable.from_db_cursor(self.cur)
        table_incidence.set_style(prettytable.PLAIN_COLUMNS)
        # print table_incidence
        table_incidence.header = True
        table_incidence.align = "l"
        table_incidence.left_padding_width = 0
        table_incidence.right_padding_width = 1
        string = table_incidence.get_string()
        string = str(string)
        result_DHF = [tuple(filter(None, map(str.strip, splitline)))
                      for line in string.splitlines()
                      for splitline in [line.split(" ")]
                      if len(splitline) > 1]
        with open('./data/incidence_data_DHF.csv', 'wb') as outcsv:
            writer = csv.writer(outcsv)
            writer.writerows(result_DHF)

        self.result_DHF = result_DHF
        self.result_DF = result_DF
        return result_DF, result_DHF
    
    def frecuency_per_day_and_week(self):
        result_DF, result_DHF = self.hemorrhagic_query()
        occurrences_dates_DF = [result_DF[i][0] for i in
                                range(1, len(result_DF))]
        occurrences_dates_DHF = [result_DHF[i][0] for i in
                                 range(1, len(result_DHF))]

        frecuency_per_date_DF = Counter(occurrences_dates_DF)
        frecuency_per_date_DHF = Counter(occurrences_dates_DHF)
        frecuency_per_date_DF_list = []
        frecuency_per_date_DHF_list = []
        # Change to list
        for key, value in frecuency_per_date_DF.iteritems():
            temp = [key, value]
            frecuency_per_date_DF_list.append(temp)

        for key, value in frecuency_per_date_DHF.iteritems():
            temp = [key, value]
            frecuency_per_date_DHF_list.append(temp)
        frecuency_per_date_week_DF_list = []
        frecuency_per_date_week_DHF_list = []
        # adding headers
        # frecuency_per_date_week_DF_list.append(['date', 'week', 'frecuency'])
        # frecuency_per_date_week_DHF_list.append(['date', 'week',
        # 'frecuency'])
        for element in frecuency_per_date_DF_list:
            dt = element[0]
            dt = datetime.datetime.strptime(dt, "%m/%d/%y").date()
            week = dt.strftime('%U')
            frecuency_per_date_week_DF_list.append([element[0], week,
                                                    element[1]])
        for element in frecuency_per_date_DHF_list:
            dt = element[0]
            dt = datetime.datetime.strptime(dt, "%m/%d/%y").date()
            week = dt.strftime('%U')
            frecuency_per_date_week_DHF_list.append([element[0], week,
                                                     element[1]])
        frecuency_per_date_week_DF_list.sort()
        frecuency_per_date_week_DHF_list.sort()
        f1 = './data/frecuency_per_day_per_week_FD.dat'
        f2 = './data/frecuency_per_day_per_week_FHD.dat'
        #
        with open(f1, 'w') as file_handler:
            file_handler.seek(0)
            file_handler.write(
                "{},{},{}\n".format('date', 'week', 'frecuency'))
            for item in frecuency_per_date_week_DF_list:
                file_handler.write(
                    "{},{},{}\n".format(item[0], item[1], item[2]))
        file_handler.close()
        #
        with open(f2, 'w') as file_handler:
            file_handler.seek(0)
            file_handler.write(
                "{},{},{}\n".format('date', 'week', 'frecuency'))
            for item in frecuency_per_date_week_DHF_list:
                file_handler.write(
                    "{},{},{}\n".format(item[0], item[1], item[2]))
        file_handler.close()

        frecuency_week_DF = [frecuency_per_date_week_DF_list[i][1]
                             for i in
                             range(len(frecuency_per_date_week_DF_list))]
        frecuency_week_DHF = [frecuency_per_date_week_DHF_list[i][1]
                              for i in
                              range(len(frecuency_per_date_week_DHF_list))]
        first_week_DF = int(frecuency_week_DF[0])
        first_week_DHF = int(frecuency_week_DHF[0])
        last_week_DF = int(frecuency_week_DF[-1])
        last_week_DHF = int(frecuency_week_DHF[-1])
        index_occurrence_week_DF = []
        index_occurrence_week_DHF = []
        #
        for i in range(first_week_DF, last_week_DF + 1):
            week = str(i).zfill(2)
            week_index = self.indices(frecuency_week_DF, week)
            index_occurrence_week_DF.append(week_index)

        for i in range(first_week_DHF, last_week_DHF + 1):
            week = str(i).zfill(2)
            week_index = self.indices(frecuency_week_DHF, week)
            index_occurrence_week_DHF.append(week_index)

        frecuency_per_week_DF = []
        frecuency_per_week_DHF = []

        week = first_week_DF
        for index in index_occurrence_week_DF:
            if len(index) != 0:
                frecuency = 0
                for k in index:
                    frecuency += frecuency_per_date_week_DF_list[k][2]
                frecuency_per_week_DF.append([week, frecuency])
            week += 1

        week = first_week_DHF
        for index in index_occurrence_week_DHF:
            if len(index) != 0:
                frecuency = 0
                for k in index:
                    frecuency += frecuency_per_date_week_DHF_list[k][2]
                frecuency_per_week_DHF.append([week, frecuency])
            week += 1
        frecuency_per_week_DF = np.array(frecuency_per_week_DF)
        frecuency_per_week_DHF = np.array(frecuency_per_week_DHF)
        self.frecuency_per_week_DF = frecuency_per_week_DF
        self.frecuency_per_week_DHF = frecuency_per_week_DHF
        np.savetxt('./data/frecuency_per_week_DF.dat',
                   frecuency_per_week_DF, fmt='%d', delimiter=',')
        np.savetxt('./data/frecuency_per_week_DHF.dat',
                   frecuency_per_week_DHF, fmt='%d', delimiter=',')

        return frecuency_per_week_DF, frecuency_per_week_DHF

    def plot_data_frecuency_per_week(self):
        frecuency_per_week_DF, frecuency_per_week_DHF = \
            self.frecuency_per_day_and_week()
        plt.plot(frecuency_per_week_DF[2: -1, 0],
                 frecuency_per_week_DF[2: -1, 1],
                 ls='--',
                 marker='o',
                 mfc='blue',
                 ms=10,
                 alpha=0.4
                 )
        plt.plot(frecuency_per_week_DHF[2: -1, 0],
                 frecuency_per_week_DHF[2: -1, 1],
                 linestyle='-.',
                 marker='o',
                 mfc='red',
                 ms=10,
                 alpha=0.4
                 )
        plt.show()
