"""Вспомогательные функции для filter_fastq: чтение, запись и проверки ридов."""


def read_fastq(input_fastq: str):
    """
    Читает fastq-файл по одному риду (4 строки за раз).

    Аргументы:
    input_fastq: str - путь до fastq-файла

    Выдаёт (yield) кортеж (name, sequence, plus_line, quality).
    yield отдаёт один рид и "ставит функцию на паузу" до следующего,
    поэтому весь файл не загружается в память.
    """
    with open(input_fastq) as file:
        while True:
            name = file.readline().strip()
            if name == '':
                # файл закончился
                break
            sequence = file.readline().strip()
            plus_line = file.readline().strip()
            quality = file.readline().strip()
            yield name, sequence, plus_line, quality


def write_read(output_file, read: tuple) -> None:
    """
    Записывает один рид (4 строки) в уже открытый файл.

    Аргументы:
    output_file: открытый на запись файл
    read: tuple - (name, sequence, plus_line, quality)

    Ничего не возвращает.
    """
    for line in read:
        output_file.write(line + '\n')


def make_bounds(bounds: tuple | float) -> tuple:
    """
    Приводит границы к виду (нижняя, верхняя).
    Если передано одно число, это верхняя граница, а нижняя равна 0.

    Аргументы:
    bounds: tuple из двух чисел или одно число

    Возвращает tuple (нижняя, верхняя).
    """
    if type(bounds) == int or type(bounds) == float:
        return (0, bounds)
    return (bounds[0], bounds[1])


def count_gc(sequence: str) -> float:
    """
    Считает GC-состав последовательности в процентах.

    Аргументы:
    sequence: str

    Возвращает float от 0 до 100 (0 для пустой последовательности).
    """
    if len(sequence) == 0:
        return 0
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence) * 100


def count_mean_quality(quality: str) -> float:
    """
    Считает среднее качество рида по шкале phred33.
    Качество одного нуклеотида = ASCII-код символа минус 33.

    Аргументы:
    quality: str - строка качества рида

    Возвращает float (0 для пустой строки).
    """
    if len(quality) == 0:
        return 0
    total = 0
    for symbol in quality:
        total += ord(symbol) - 33
    return total / len(quality)


def is_gc_ok(sequence: str, gc_bounds: tuple) -> bool:
    """
    Проверяет, что GC-состав попадает в gc_bounds (границы включаются).

    Аргументы:
    sequence: str
    gc_bounds: tuple - (нижняя, верхняя) в процентах

    Возвращает bool.
    """
    gc = count_gc(sequence)
    return gc_bounds[0] <= gc <= gc_bounds[1]


def is_length_ok(sequence: str, length_bounds: tuple) -> bool:
    """
    Проверяет, что длина последовательности попадает в length_bounds
    (границы включаются).

    Аргументы:
    sequence: str
    length_bounds: tuple - (нижняя, верхняя)

    Возвращает bool.
    """
    length = len(sequence)
    return length_bounds[0] <= length <= length_bounds[1]


def is_quality_ok(quality: str, quality_threshold: float) -> bool:
    """
    Проверяет, что среднее качество рида не ниже quality_threshold.

    Аргументы:
    quality: str - строка качества рида
    quality_threshold: float - порог (phred33)

    Возвращает bool.
    """
    return count_mean_quality(quality) >= quality_threshold
