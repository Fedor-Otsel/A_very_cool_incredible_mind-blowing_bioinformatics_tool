"""Вспомогательные функции для чтения fasta-файлов и результатов BLAST."""


def read_fasta(input_fasta: str):
    """
    Читает fasta-файл, в котором последовательность может быть
    разбита на несколько строк, и склеивает эти строки.

    Аргументы:
    input_fasta: str - путь до fasta-файла

    Выдаёт (yield) кортеж (header, sequence) по одной записи за раз.
    """
    header = ''
    sequence = ''
    with open(input_fasta) as file:
        for line in file:
            line = line.strip()
            if line == '':
                continue
            if line.startswith('>'):
                # началась новая запись - отдаём предыдущую
                if header != '':
                    yield header, sequence
                header = line
                sequence = ''
            else:
                sequence += line
    # последняя запись в файле
    if header != '':
        yield header, sequence


def make_oneline_name(input_fasta: str) -> str:
    """
    Придумывает имя выходного файла по имени входного.
    Пример: 'seqs.fasta' -> 'seqs_oneline.fasta'

    Аргументы:
    input_fasta: str

    Возвращает str.
    """
    if input_fasta.endswith('.fasta'):
        return input_fasta[:-6] + '_oneline.fasta'
    return input_fasta + '_oneline.fasta'


def read_blast_descriptions(input_file: str) -> list:
    """
    Читает txt-файл BLAST и для каждого запроса берёт
    первую строку из столбца Description (лучшее совпадение).

    Аргументы:
    input_file: str - путь до результатов BLAST (txt)

    Возвращает list со строками описаний.
    """
    with open(input_file) as file:
        lines = file.readlines()

    descriptions = []
    for i in range(len(lines)):
        if lines[i].startswith('Sequences producing significant alignments'):
            # lines[i + 1] - строка заголовка со словом 'Scientific'
            # lines[i + 2] - строка заголовка со словом 'Description'
            # lines[i + 3] - первое (лучшее) совпадение
            # Столбец Description заканчивается там, где начинается 'Scientific'
            column_end = lines[i + 1].index('Scientific')
            best_hit = lines[i + 3][:column_end].strip()
            descriptions.append(best_hit)
    return descriptions


def write_lines(output_file: str, lines: list) -> None:
    """
    Записывает строки в файл, каждую с новой строки.

    Аргументы:
    output_file: str - имя выходного файла
    lines: list - список строк

    Ничего не возвращает.
    """
    with open(output_file, 'w') as file:
        for line in lines:
            file.write(line + '\n')
