"""
Здесь много дублирования кода, что может вызвать страдания у любителей программирования
Увы в этом курсе плевать на чистоту кода
"""

from copy import deepcopy

from typing import Optional, Tuple

DEBUG = False
EPS = 1e-6

DATAPATH = r"D:\develop\IU7-CA\lab_01\data\testset.txt"
DATA1PATH = r"D:\develop\IU7-CA\lab_01\data\func1.txt"
DATA2PATH = r"D:\develop\IU7-CA\lab_01\data\func2.txt"


class DataSlice:
    def __init__(self, x: float, y: float,
                 derivative: Optional[float] = None,
                 double_derivative: Optional[float] = None) -> None:
        self.x: float = x
        self.y: float = y
        self.derivative: float | None = derivative
        self.double_derivative: float | None = double_derivative


class Model:
    def __init__(self) -> None:
        self.data_list: list[DataSlice] = []
        self.newton_deg: int = 5
        self.ermite_nodes: int = 5

        self.data1: list[DataSlice] = []
        self.data2: list[DataSlice] = []

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
        while pos < len(self.data_list) and self.data_list[pos].x < x:
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

        coeff_table: list[list[float]] = [
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
            xcoef *= x - coeff_table[0][i]
            res += coeff_table[i+2][0] * xcoef

        return res

    def approx_by_ermite(self, x):
        if self.ermite_nodes == 0:
            return print("Не определено кол-во точек для полинома Эрмита")
        if len(self.data_list) < self.ermite_nodes:
            return print("Для данного количества точек полинома Эрмита недостаточно входных данных функции")

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

        nodes: list[DataSlice] = self.data_list[index_low: index_high]

        coeffs = self.get_ermite_coeffs(nodes)

        n = len(coeffs[0])

        res = coeffs[1][0]
        xcoef: float = 1.0

        for i in range(n-1):
            xcoef *= x - coeffs[0][i]
            res += coeffs[i+2][0] * xcoef

        return res

    def get_ermite_coeffs(self, nodes: list[DataSlice]) -> list[float]:
        n: int = 0
        for slice in nodes:
            n += 1
            if slice.derivative is not None:
                n += 1
            if slice.double_derivative is not None:
                n += 1

        table: list[list[float | None]] = \
            [
                [None for _ in range(n)],
                [None for _ in range(n)]
        ]

        for i in range(n - 1):
            table.append([None for _ in range(n - i - 1)])

        cur: int = 0
        for slice in nodes:
            table[0][cur] = slice.x
            table[1][cur] = slice.y
            cur += 1
            if slice.derivative is not None:
                table[0][cur] = slice.x
                table[1][cur] = slice.y
                table[2][cur - 1] = slice.derivative
                cur += 1
            if slice.double_derivative is not None:
                table[0][cur] = slice.x
                table[1][cur] = slice.y
                table[2][cur - 1] = slice.derivative
                table[3][cur - 2] = slice.double_derivative
                cur += 1

        for i in range(n - 1):
            for j in range(n - i - 1):
                if table[i + 2][j] is None:
                    table[i + 2][j] = (
                        (table[i+1][j] - table[i+1][j+1]) /
                        (table[0][j] - table[0][j + i + 1])
                    )

        if DEBUG:
            for line in table:
                for item in line:
                    if item is None:
                        print("None", end="\t")
                    else:
                        print(f"{item:0.4}", end="\t")
                print("")

        return table

    def get_root_by_newton(self):
        if self.newton_deg == 0:
            return print("Не определена степень полинома Ньютона")
        if len(self.data_list) < self.newton_deg+1:
            return print("Для данной степени полинома ньютона недостаточно входных данных функции")

        x = 0

        pos = 0
        minmod = abs(self.data_list[pos].y)

        for i in range(len(self.data_list)):
            if abs(self.data_list[i].y) < minmod:
                minmod = abs(self.data_list[i].y)
                pos = i

        index_low = pos - (self.newton_deg + 1) // 2
        index_high = pos + (self.newton_deg + 1) // 2 + \
            (self.newton_deg + 1) % 2

        if index_low < 0:
            index_high -= index_low
            index_low = 0

        if index_high >= len(self.data_list) - 1:
            index_low -= index_high - len(self.data_list)
            index_high = len(self.data_list)

        coeff_table: list[list[float]] = [
            [s.y for s in self.data_list[index_low: index_high]],
            [s.x for s in self.data_list[index_low: index_high]],
        ]

        for _ in range(self.newton_deg):
            coeff_table.append([])

        for i in range(self.newton_deg):
            for j in range(self.newton_deg - i):
                coeff_table[i+2].append(
                    (coeff_table[i+1][j] - coeff_table[i+1][j+1]) /
                    (coeff_table[0][j] - coeff_table[0][j+i+1])
                )

        if True:
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
            xcoef *= x - coeff_table[0][i]
            res += coeff_table[i+2][0] * xcoef

        return res

    def get_root_by_ermite(self):
        if self.ermite_nodes == 0:
            return print("Не определено кол-во точек для полинома Эрмита")
        if len(self.data_list) < self.ermite_nodes:
            return print("Для данного количества точек полинома Эрмита недостаточно входных данных функции")

        x = 0

        pos = 0
        minmod = abs(self.data_list[pos].x)

        for i in range(len(self.data_list)):
            if abs(self.data_list[pos].x) < minmod:
                minmod = abs(self.data_list[pos].x)
                pos = i

        index_low = pos - (self.ermite_nodes) // 2
        index_high = pos + (self.ermite_nodes) // 2 + \
            (self.ermite_nodes) % 2

        if index_low < 0:
            index_high -= index_low
            index_low = 0

        if index_high >= len(self.data_list) - 1:
            index_low -= index_high - len(self.data_list) + 1
            index_high = len(self.data_list) - 1

        nodes: list[DataSlice] = deepcopy(
            self.data_list[index_low: index_high])

        for node in nodes:
            tmp = node.x
            node.x = node.y
            node.y = tmp

            if node.derivative is not None:
                node.derivative = 1 / \
                    node.derivative if abs(
                        node.derivative) > EPS else None

            if node.double_derivative is not None:
                node.double_derivative = 1 / \
                    node.double_derivative if abs(
                        node.double_derivative) > EPS else None

        coeffs = self.get_ermite_coeffs(nodes)

        n = len(coeffs[0])

        res = coeffs[1][0]
        xcoef: float = 1.0

        for i in range(n-1):
            xcoef *= x - coeffs[0][i]
            res += coeffs[i+2][0] * xcoef

        return res

    def find_system_root(self) -> Tuple[float, float]:
        res: list[DataSlice] = []

        for i in [i/10 for i in range(1, 11)]:
            x = i
            self.data_list = self.data2
            y = self.approx_by_newton(x)
            self.data_list = self.data1
            y -= self.approx_by_newton(x)

            res.append(DataSlice(x, y, None, None))

        for slice in res:
            print(slice.x, slice.y, slice.derivative,
                  slice.double_derivative, sep="\t")

        self.data_list = res
        x = self.get_root_by_newton()
        self.data_list = self.data1

        return (x, self.approx_by_newton(x))


class Controller:
    def __init__(self, model: Model) -> None:
        self.model: Model = model


class View:
    def __init__(self, model: Model, controller: Controller) -> None:
        self.model: Model = model
        self.controller: Controller = controller
        self.running = False
        self.options = {
            "1": self.read_from_file,       # read data
            "2": self.read_newton_deg,      # newton deg
            "3": self.read_ermite_nodes,    # ermite nodes
            "4": self.approx_by_newton,     # newton by x
            "5": self.approx_by_ermite,     # ermite by x
            "6": self.get_roots,            # roots
            "7": self.test_both,                       # table
            "8": self.get_sys_root,
        }

    def get_roots(self):
        newton_root: float | None = self.model.get_root_by_newton()
        ermite_root: float | None = self.model.get_root_by_ermite()

        if newton_root is not None:
            print(f"Корень по Ньютону: {newton_root:.6}")
        if ermite_root is not None:
            print(f"Корень по Эрмиту:  {ermite_root:.6}")

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
        filename = DATAPATH if DATAPATH != "" else read_filename()
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

    def read_funcs(self):
        try:
            with open(DATA1PATH, "r") as f:
                lines = f.readlines()
                for line in lines:
                    l = parse_line(line)
                    if l is None:
                        continue
                    self.model.data1.append(l)
        except OSError:
            return print("Не удалось открыть файл")

        try:
            with open(DATA2PATH, "r") as f:
                lines = f.readlines()
                for line in lines:
                    l = parse_line(line)
                    if l is None:
                        continue
                    self.model.data2.append(l)
        except OSError:
            return print("Не удалось открыть файл")

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
        if self.model.ermite_nodes > len(self.model.data_list):
            return print("Кол-во узлов больше количества входных данных")

        try:
            x = float(input("Введите значение x: "))
        except ValueError:
            return print("Не удалось считать значение x")

        res = self.model.approx_by_ermite(x)
        print(f"Результат аппроксимации: {res:.6}")

        return res

    def test_both(self):
        if not len(self.model.data_list):
            return print("Сначала добавьте входные данные")
        if 5+1 > len(self.model.data_list):
            return print("Степень полинома Ньютона больше количества входных данных")
        if 5 > len(self.model.data_list):
            return print("Кол-во узлов больше количества входных данных")

        try:
            x = float(input("Введите значение x: "))
        except ValueError:
            return print("Не удалось считать значение x")

        print(f"{'N':^3}|{'Ньютон':^20}|{'Эрмит':^20}")

        for i in range(5):
            self.model.ermite_nodes = i + 1
            self.model.newton_deg = i + 1

            new = self.model.approx_by_newton(x)
            erm = self.model.approx_by_ermite(x)

            print(f"""{i+1:^3}|{f"{new:.6}":^20}|{f"{erm:.6}":^20}""")

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
                  "8: Найти корень системы\n" +
                  "0: Выход")
            option = input("> ")

            if option == "0":
                self.running = False
            elif len(option) > 1 or option not in "12345678":
                continue
            else:
                self.options[option]()

    def run(self):
        self.running = True
        self.menu_loop()

    def get_sys_root(self):
        self.read_funcs()

        x, y = self.model.find_system_root()

        print(x, y)


def read_filename() -> str:
    return input("Введите имя файла: ")


def parse_line(line: str) -> DataSlice | None:
    data = [*map(float, line.strip().split())]
    if len(data) < 2 or len(data) > 4:
        if DEBUG:
            print("Строка пропущена: ", *data)
        return None

    return DataSlice(*data)


if __name__ == "__main__":
    model: Model = Model()
    controller: Controller = Controller(model)
    view: View = View(model, controller)

    if (DATAPATH):
        view.read_from_file()

    view.run()
