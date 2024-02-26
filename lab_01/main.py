from typing import *
from copy import copy
import numpy as np

DEBUG = True


class DataSlice:
    def __init__(self, x: float, y: float,
                 derivative: float,
                 double_derivative: float) -> None:
        self.x: float = x
        self.y: float = y
        self.derivative: float = derivative
        self.double_derivative: float = double_derivative


class Model:
    def __init__(self) -> None:
        self.data_list: list[DataSlice] = []
        self.newton_deg: int = 0
        self.ermite_nodes: int = 0

    def push_slice(self, data: DataSlice):
        if data:
            self.data_list.append(data)
            if DEBUG:
                print("добавлено: ", data.x, data.y,
                      data.derivative, data.double_derivative)

    def sort_slices(self):
        self.data_list.sort(key=lambda s: s.x)

    def approx_by_newton(self, x):
        if self.newton_deg == 0:
            return print("Не определена степень полинома Ньютона")
        if len(self.data_list) < self.newton_deg+1:
            return print("Для данной степени полинома ньютона недостаточно входных данных функции")

        pos = 0
        while self.data_list[pos].x < x and pos < len(self.data_list):
            pos += 1

        index_low = pos - (self.newton_deg + 1) // 2
        index_high = pos + (self.newton_deg + 1) // 2 + \
            (self.newton_deg + 1) % 2

        if index_low < 0:
            index_high -= index_low
            index_low = 0

        if index_high >= len(self.data_list) - 1:
            index_low -= index_high - len(self.data_list) + 1
            index_high = len(self.data_list) - 1

        coeff_table: List[List[float]] = [
            [s.x for s in self.data_list[index_low: index_high]],
            [s.y for s in self.data_list[index_low: index_high]],
        ]

        for _ in range(self.newton_deg):
            coeff_table.append([])

        for i in range(self.newton_deg):
            for j in range(self.newton_deg - i):
                coeff_table[i+2].append(
                    (coeff_table[i+1][j] - coeff_table[i+1][j+1]) /
                    (coeff_table[0][j] - coeff_table[0][j+i+1])
                )

        if DEBUG:
            print("Таблица коэф-тов:")
            for i, line in enumerate(coeff_table):
                if (i < 2):
                    print(" x:" if i == 0 else " y:", end='\t')
                else:
                    print(f"y{i-2}:", end='\t')

                for item in line:
                    print(f"{item:.3}", end='\t')
                print()

        res = coeff_table[1][0]
        xcoef: float = 1.0

        for i in range(self.newton_deg):
            # if DEBUG:
            #     print(f"xcoef {i} : {xcoef:.3} * ({x} - {coeff_table[0][i]})")
            xcoef *= x - coeff_table[0][i]
            res += coeff_table[i+2][0] * xcoef

        return res

    def transpose_matrix(self, matrix):
        return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]

    def matrix_minor(self, matrix, i, j):
        return [row[:j] + row[j+1:] for row in (matrix[:i]+matrix[i+1:])]

    def matrix_determinant(self, matrix):
        if len(matrix) == 2:
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]

        determinant = 0
        for c in range(len(matrix)):
            determinant += ((-1) ** c) * matrix[0][c] * \
                self.matrix_determinant(self.matrix_minor(matrix, 0, c))
        return determinant

    def matrix_inverse(self, matrix):
        determinant = self.matrix_determinant(matrix)
        if len(matrix) == 2:
            return [[matrix[1][1] / determinant, -1 * matrix[0][1] / determinant],
                    [-1 * matrix[1][0] / determinant, matrix[0][0] / determinant]]

        cofactors = []
        for r in range(len(matrix)):
            cofactor_row = []
            for c in range(len(matrix)):
                minor = self.matrix_minor(matrix, r, c)
                cofactor_row.append(((-1) ** (r + c)) *
                                    self.matrix_determinant(minor))
            cofactors.append(cofactor_row)
        cofactors = self.transpose_matrix(cofactors)
        for r in range(len(cofactors)):
            for c in range(len(cofactors)):
                cofactors[r][c] = cofactors[r][c] / determinant
        return cofactors

    def approx_by_ermite(self, x):
        pos = 0
        while self.data_list[pos].x < x and pos < len(self.data_list):
            pos += 1

        index_low = pos - (self.ermite_nodes) // 2
        index_high = pos + (self.ermite_nodes) // 2 + \
            (self.ermite_nodes) % 2

        if index_low < 0:
            index_high -= index_low
            index_low = 0

        if index_high >= len(self.data_list) - 1:
            index_low -= index_high - len(self.data_list) + 1
            index_high = len(self.data_list) - 1

        coeff_table: List[List[float]] = [
            [s.x for s in self.data_list[index_low: index_high]],
            [s.y for s in self.data_list[index_low: index_high]],
            [s.derivative for s in self.data_list[index_low: index_high]],
            [s.double_derivative for s in self.data_list[index_low: index_high]],
        ]

        coeffs = self.get_ermite_coeffs(coeff_table)

        sum = 0

        print(coeffs)

        for i in range(len(coeffs)):
            sum += coeffs[i] * (x ** (len(coeffs) - i - 1))

        return sum

    def get_ermite_coeffs(self, coeff_table: list[list[float]]) -> list[float]:
        n = len(coeff_table[0])*3
        xlen = len(coeff_table[0])

        matrix = [[0.0 for _ in range(n)] for _ in range(n)]
        b_col = [[0.0] for _ in range(n)]

        for i in range(xlen):
            for j in range(n):
                if (n - j - 1) < 0:
                    continue
                matrix[i][j] = coeff_table[0][i] ** (n - j - 1)
            b_col[i][0] = coeff_table[1][i]

        for i in range(xlen):
            for j in range(n):
                if (n - j - 2) < 0:
                    continue
                matrix[i +
                       xlen][j] = (n - j - 1) * coeff_table[0][i] ** (n - j - 2)
            b_col[i + xlen][0] = coeff_table[2][i]

        for i in range(xlen):
            for j in range(n):
                if (n - j - 3) < 0:
                    continue
                matrix[i +
                       xlen * 2][j] = (n - j - 1) * (n - j - 2) * coeff_table[0][i] ** (n - j - 3)
            b_col[i + xlen * 2][0] = coeff_table[3][i]

        tmp = np.matrix(matrix)

        inv = np.linalg.inv(tmp)
        b = np.matrix(b_col)

        result = np.matmul(inv, b)
        result = result.tolist()

        return [result[i][0] for i in range(n)]

    def matrix_mul(self, m1: list[list[float]], m2: list[list[float]]) -> list[list[float]]:
        m = len(m1)                                            # a: m × n
        n = len(m2)                                            # b: n × k
        k = len(m2[0])

        res = [[0.0 for _ in range(k)] for _ in range(m)]    # c: m × k

        for i in range(m):
            for j in range(k):
                res[i][j] = sum(m1[i][kk] * m2[kk][j] for kk in range(n))

        return res


class Controller:
    def __init__(self, model: Model) -> None:
        self.model: Model = model


class View:
    def __init__(self, model: Model, controller: Controller) -> None:
        self.model: Model = model
        self.controller: Controller = controller
        self.running = False
        self.options = {
            "1": self.read_from_file,  # read data
            "2": self.read_newton_deg,  # newton deg
            "3": self.read_ermite_nodes,  # ermite nodes
            "4": self.approx_by_newton,  # newton by x
            "5": self.approx_by_ermite,  # ermite by x
            "6": ...,  # roots
            "7": ...,  # table
        }

    def print_status(self):
        print("Текущий статус:\n" +
              f"Количество строк в таблице: {len(self.model.data_list)}\n" +
              f"Степень полинома Ньютона: {self.model.newton_deg if self.model.newton_deg else 'Не задано'}\n" +
              f"Кол-во узлов полинома Эрмита: {self.model.ermite_nodes if self.model.ermite_nodes else 'Не задано'}"
              )

    def read_newton_deg(self):
        try:
            newdeg = int(input("Введите степень полинома Ньютона: "))
        except ValueError:
            return print("Ошибка считывания")

        if newdeg <= 0:
            return print("Степень полинома слишком маленькая")

        if newdeg >= 7:
            print("Степень полинома больше 7 (предупреждение)")

        self.model.newton_deg = newdeg

    def read_ermite_nodes(self):
        try:
            newnodes = int(input("Введите Количество узлов полинома Эрмита: "))
        except ValueError:
            return print("Ошибка считывания")

        if newnodes <= 0:
            return print("Кол-во узлов слишком мало")

        if newnodes >= 7:
            print("Кол-во узлов больше 7 (предупреждение)")

        self.model.ermite_nodes = newnodes

    def read_from_file(self):
        filename = read_filename()
        try:
            with open(filename, "r") as f:
                lines = f.readlines()
                for line in lines:
                    l = parse_line(line)
                    if l is None:
                        continue
                    self.model.push_slice(l)
        except OSError:
            return print("Не удалось открыть файл")

        if len(self.model.data_list) == 0:
            return print("Не удалось считать ни одной строки с данными")

        self.model.sort_slices()

    def approx_by_newton(self):
        if not self.model.newton_deg:
            return print("Сначала введите степень полинома")
        if not len(self.model.data_list):
            return print("Сначала добавьте входные данные")
        if self.model.newton_deg+1 > len(self.model.data_list):
            return print("Степень полинома Ньютона больше количества входных данных")

        try:
            x = float(input("Введите значение x: "))
        except ValueError:
            return print("Не удалось считать значение x")

        res = self.model.approx_by_newton(x)
        print(f"Результат аппроксимации: {res:.6}")

        return res

    def approx_by_ermite(self):
        if not self.model.ermite_nodes:
            return print("Сначала введите кол-во узлов")
        if not len(self.model.data_list):
            return print("Сначала добавьте входные данные")
        if self.model.newton_deg+1 > len(self.model.data_list):
            return print("Кол-во узлов больше количества входных данных")

        try:
            x = float(input("Введите значение x: "))
        except ValueError:
            return print("Не удалось считать значение x")

        res = self.model.approx_by_ermite(x)
        print(f"Результат аппроксимации: {res:.6}")

        return res

    def menu_loop(self):
        while self.running:
            self.print_status()
            print("Меню:\n" +
                  "1: Считать данные из файла\n" +
                  "2: Задать степень полинома Ньютона\n" +
                  "3: Задать количество узлов полинома Эрмита\n" +
                  "4: Получить приблизительное значение y(x) с помощью полинома Ньютона\n" +
                  "5: Получить приблизительное значение y(x) с помощью полинома Эрмита\n" +
                  "6: Найти корень с помощью обоих полиномов\n" +
                  "7: Построить таблицу сранения полиномов при заданном x\n" +
                  "0: Выход")
            option = input("> ")

            if option == "0":
                self.running = False
            elif len(option) > 1 or option not in "1234567":
                continue
            else:
                self.options[option]()

    def run(self):
        self.running = True
        self.menu_loop()


def read_filename() -> str:
    return input("Введите имя файла: ")


def parse_line(line: str) -> DataSlice | None:
    data = [*map(float, line.strip().split())]
    if len(data) != 4:
        if DEBUG:
            print("Строка пропущена: ", *data)
        return None

    return DataSlice(*data)


if __name__ == "__main__":
    model: Model = Model()
    controller: Controller = Controller(model)
    view: View = View(model, controller)

    view.run()
