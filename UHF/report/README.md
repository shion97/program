# 実験01〜08のレポート

各 `experiment1/` 〜 `experiment8/` に `report.tex` とコンパイル済み `report.pdf` を置く。
原資料は `../results_vault/experiment_01/` 〜 `experiment_08/` と `../experiment_conditions.md`。
計算結果は変更せず、各実験内の全PNG（再利用ケースを含む計232枚）を参照する。

章構成は生成物に対応する7章：計算ログ、密度CSV、SCFエネルギー、SCF残差、軌道エネルギー、局所磁化、境界条件によるエネルギー差。
記録された最終値と密度CSVを根拠に、条件依存性とOBC/PBC差を説明する。
未収束値は暫定値と明記し、境界エネルギー差の確定比較から除外する。
軌道固有値の数値ファイルはないため、保存図を用い、未収束密度から固有値を再計算して保存図と同一視しない。

## 完了確認

- 8本のLaTeXがコンパイルでき、日本語と図が表示される。
- 原資料の全232枚が各レポートに1回ずつ含まれ、参照切れがない。
- 全56ケースの最終数値、収束状態、粒子数を原資料と照合する。
- 各実験で全体傾向、変化させた条件の影響、境界条件の違いと解釈上の限界を説明する。

## 再コンパイル

各experimentディレクトリで `latexmk -norc -lualatex -interaction=nonstopmode -halt-on-error report.tex` を実行する。
TeX LiveのLuaLaTeX、luatexja、Harano Ajiフォントを使用する。
画像は元のresults_vaultを相対パスで参照するため、reportとresults_vaultの位置関係を保つ。
補助ファイルと検証記録は `_work/` に置く。

## 成果物

| 実験 | PDF | LaTeX |
|---|---|---|
| 1 | [report.pdf](experiment1/report.pdf) | [report.tex](experiment1/report.tex) |
| 2 | [report.pdf](experiment2/report.pdf) | [report.tex](experiment2/report.tex) |
| 3 | [report.pdf](experiment3/report.pdf) | [report.tex](experiment3/report.tex) |
| 4 | [report.pdf](experiment4/report.pdf) | [report.tex](experiment4/report.tex) |
| 5 | [report.pdf](experiment5/report.pdf) | [report.tex](experiment5/report.tex) |
| 6 | [report.pdf](experiment6/report.pdf) | [report.tex](experiment6/report.tex) |
| 7 | [report.pdf](experiment7/report.pdf) | [report.tex](experiment7/report.tex) |
| 8 | [report.pdf](experiment8/report.pdf) | [report.tex](experiment8/report.tex) |

`numerical_summary.json` は各レポートの数値表の根拠となる集計値。保存図とログの数値は原資料を参照する。

制限された実行環境でLuaLaTeXのキャッシュが書き込めない場合は、書き込み可能な `report/_work/texcache` を作成し、当該プロセスの `TEXMFVAR` と `TEXMFCACHE` にその絶対パスを設定する。
