"""Bio File Tools - простые утилиты для fastq, fasta и результатов BLAST."""

import os

from helpers.fastq import read_fastq, write_read, make_bounds
from helpers.fastq import is_gc_ok, is_length_ok, is_quality_ok
from helpers.bio_files import read_fasta, make_oneline_name
from helpers.bio_files import read_blast_descriptions, write_lines


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple | float = (0, 100),
    length_bounds: tuple | int = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """
    Фильтрует риды из fastq-файла по GC-составу, длине и среднему качеству.
    Риды проверяются по одному, и подходящие сразу записываются
    в файл filtered/<output_fastq>.

    Аргументы:
    input_fastq: str - путь до входного fastq-файла
    output_fastq: str - имя выходного файла (сохраняется в папку filtered)
    gc_bounds: tuple (нижняя, верхняя) или одно число (верхняя граница),
        в процентах. По умолчанию (0, 100)
    length_bounds: tuple (нижняя, верхняя) или одно число (верхняя граница).
        По умолчанию (0, 2**32)
    quality_threshold: float - минимальное среднее качество (phred33).
        По умолчанию 0

    Все границы включаются.
    Ничего не возвращает.
    """
    gc_bounds = make_bounds(gc_bounds)
    length_bounds = make_bounds(length_bounds)

    # создаём папку filtered, если её ещё нет
    os.makedirs('filtered', exist_ok=True)
    output_path = os.path.join('filtered', output_fastq)

    with open(output_path, 'w') as output_file:
        for read in read_fastq(input_fastq):
            sequence = read[1]
            quality = read[3]
            if not is_gc_ok(sequence, gc_bounds):
                continue
            if not is_length_ok(sequence, length_bounds):
                continue
            if not is_quality_ok(quality, quality_threshold):
                continue
            write_read(output_file, read)


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str = None
) -> None:
    """
    Записывает каждую последовательность fasta-файла в одну строку.

    Аргументы:
    input_fasta: str - путь до входного fasta-файла
    output_fasta: str - имя выходного файла (необязательно).
        Если не указано, файл называется '<имя входного>_oneline.fasta'

    Ничего не возвращает.
    """
    if output_fasta is None:
        output_fasta = make_oneline_name(input_fasta)

    with open(output_fasta, 'w') as output_file:
        for header, sequence in read_fasta(input_fasta):
            output_file.write(header + '\n')
            output_file.write(sequence + '\n')


def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Для каждого запроса в txt-файле BLAST берёт лучшее совпадение
    (первую строку столбца Description) и сохраняет список белков
    одним столбцом, без повторов, по алфавиту.

    Аргументы:
    input_file: str - путь до результатов BLAST (txt)
    output_file: str - имя выходного файла

    Ничего не возвращает.
    """
    descriptions = read_blast_descriptions(input_file)
    # set убирает повторы, key=str.lower - сортировка без учёта регистра
    proteins = sorted(set(descriptions), key=str.lower)
    write_lines(output_file, proteins)
