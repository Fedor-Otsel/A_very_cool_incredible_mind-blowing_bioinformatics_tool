# Bio File Tools

Небольшой пакет на Python для работы с биоинформатическими файлами: **fastq**, **fasta** и результатами **BLAST**.
Работает только на стандартной библиотеке, ничего ставить не нужно.

## Установка

Нужен Python 3.10 или новее.

```bash
git clone https://github.com/<your-username>/bio-file-tools.git
cd bio-file-tools
```

Функции импортируются из главного скрипта `bio_file_tools.py` (запускайте Python из папки репозитория):

```python
from bio_file_tools import filter_fastq, convert_multiline_fasta_to_oneline, parse_blast_output
```

## Что умеет пакет

| Функция | Что делает |
|---|---|
| `filter_fastq` | Фильтрует риды из fastq по GC-составу, длине и среднему качеству |
| `convert_multiline_fasta_to_oneline` | Склеивает многострочные последовательности fasta в одну строку |
| `parse_blast_output` | Достаёт лучшее совпадение для каждого запроса из txt-выдачи BLAST |

### `filter_fastq`

```python
filter_fastq(input_fastq, output_fastq, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0)
```

- `input_fastq` — путь до входного fastq-файла
- `output_fastq` — имя выходного файла; результат сохраняется в папку `filtered/` (папка создаётся сама)
- `gc_bounds` — границы GC-состава в процентах. Одно число = верхняя граница
- `length_bounds` — границы длины рида. Одно число = верхняя граница
- `quality_threshold` — минимальное среднее качество рида (phred33)

Обе границы интервалов включаются. Файл читается по одному риду: прочитали → проверили → сразу записали, поэтому в памяти не хранится весь файл.

```python
filter_fastq('example_fastq.fastq', 'good_reads.fastq',
             gc_bounds=(20, 80), length_bounds=(50, 100), quality_threshold=35)
# -> filtered/good_reads.fastq

filter_fastq('example_fastq.fastq', 'low_gc.fastq', gc_bounds=44.4)
# риды с GC-составом не больше 44.4%
```

### `convert_multiline_fasta_to_oneline`

```python
convert_multiline_fasta_to_oneline(input_fasta, output_fasta=None)
```

Читает fasta-файл, где последовательность разбита на несколько строк, и сохраняет новый файл, где каждая последовательность записана в одну строку. Если `output_fasta` не указан, файл называется `<имя входного>_oneline.fasta`.

```python
convert_multiline_fasta_to_oneline('example_multiline_fasta.fasta')
# -> example_multiline_fasta_oneline.fasta
```

### `parse_blast_output`

```python
parse_blast_output(input_file, output_file)
```

Читает txt-файл с результатами BLAST. Для каждого запроса (блок `Sequences producing significant alignments:`) берёт первую строку столбца `Description`, то есть лучшее совпадение. Белки записываются в один столбец, без повторов, по алфавиту (без учёта регистра). Описания сохраняются так, как они записаны в таблице BLAST, поэтому длинные названия могут заканчиваться на `...`.

```python
parse_blast_output('example_blast_results.txt', 'best_hits.txt')
```

## Структура

```
bio-file-tools/
|- README.md
|- bio_file_tools.py   # главный скрипт: filter_fastq и функции части B
|- helpers/
    |- fastq.py        # чтение/запись fastq и проверки ридов
    |- bio_files.py    # чтение fasta и BLAST, запись результатов
```

## Автор

Fedor — по вопросам и предложениям пишите в Issues репозитория.
