from typing import *
from copy import copy

DEBUG = True


class DataSlice:
    def __init__(self, x: float, y: float,
                 derivative: Optional[float] = None,
                 double_derivative: Optional[float] = None) -> None:
        self.x: float = x
        self.y: float = y
        self.derivative: float = derivative
        self.double_derivative: float = double_derivative


class Model:
    def __init__(self) -> None:
        self.data_list: list[DataSlice] = []
        self.newton_deg: int = None
        self.ermite_nodes: int = None

    def push_slice(self, data: DataSlice):
        if data:
            self.data_list.append(data)
            if DEBUG:
                print("добавлено: ", data.x, data.y,
                      data.derivative, data.double_derivative)

    def sort_slices(self):
        self.data_list.sort(key=lambda s: s.x)

    def approx_by_newton(self, x):
        if self.newton_deg is None:
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
            "3": ...,  # ermite nodes
            "4": self.approx_by_newton,  # newton by x
            "5": ...,  # ermite by x
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
            print("Степень полинома больше 7, это может привести к неточностям")

        self.model.newton_deg = newdeg

    def read_from_file(self):
        filename = read_filename()
        try:
            with open(filename, "r") as f:
                lines = f.readlines()
                for line in lines:
                    self.model.push_slice(parse_line(line))
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


def parse_line(line: str) -> DataSlice:
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
