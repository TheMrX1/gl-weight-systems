# Конвей и Киношита–Терасака

Для двух слов из [results.txt](results.txt) получены полные универсальные gl-веса над целыми числами. В каждом многочлене 587 ненулевых мономов. Коэффициенты совпадают, [разность равна нулю](difference.txt). Второе слово получено мутацией диаграммы Конвея из [статьи Чмутова и Ландо](https://arxiv.org/abs/0704.1313).

| Диаграмма | Многочлен | Коэффициенты для программ | Вход |
|---|---|---|---|
| Конвей | [TXT](conway.txt) | [JSON](conway.json) | [Хордовое слово](conway.word) |
| Киношита–Терасака | [TXT](kinoshita-terasaka.txt) | [JSON](kinoshita-terasaka.json) | [Хордовое слово](kinoshita-terasaka.word) |

В JSON массив степеней упорядочен как N, C1, …, C22. Коэффициент хранится десятичной строкой.

Исходные остатки и значения Казимиров записаны в файлах с суффиксом '_numeric.json'. Файлы с 'hook' и 'square' содержат исходные данные свёрток для столбцов 0 и 7.

## Запуск

Вычислитель находится в [cpp/src/symbolic_gl_low_memory.cpp](../../cpp/src/symbolic_gl_low_memory.cpp). После конфигурации CMake по инструкции в корне репозитория из той же папки выполните:

~~~sh
set -e
cmake --build build --target symbolic_low_memory
out=build/gl-pair
mkdir -p "$out"
for name in conway kinoshita-terasaka; do
  word="results/kinoshita-terasaka-and-conway/$name.word"
  checkpoint="$out/$name.dag"
  CHORD_DAG_CHECKPOINT="$checkpoint" CHORD_DAG_BUILD_ONLY=1 \
    build/symbolic_low_memory --word < "$word" &&
  CHORD_DAG_CHECKPOINT="$checkpoint" \
    build/symbolic_low_memory --word < "$word" > "$out/$name.json.partial"
  mv "$out/$name.json.partial" "$out/$name.json"
done
cmp "$out/conway.json" "$out/kinoshita-terasaka.json"
python3 scripts/verify_conway_kinoshita.py
~~~

Последняя команда проверяет сохранённые коэффициенты и подстановки. Временные графы для воспроизведения создаются в папке сборки.
