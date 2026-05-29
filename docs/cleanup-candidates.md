# Cleanup Candidates

このリストは、コードレビュー時に見つけた「不要と思われる」「未使用」「到達不能」「冗長」コードの候補です。
**いずれも自動修正は行っていません。** 削除・修正の判断はユーザーの確認後に行ってください。

判定の確信度: H = 高 / M = 中 / L = 低（呼び出し元が見つからないだけで、将来使う可能性あり）

---

## 1. 完全に未使用のモジュール / シンボル

### 1.1 `src/acmg_classifier/io/json_writer.py` — モジュール全体 (H)
- `write_json()` がプロジェクト内のどこからも import されていない。
- TSV 出力 (`tsv_writer.py`) のみが pipeline / cli から使われている。
- JSON 出力機能を将来提供する予定がなければ、モジュールごと削除可能。

### 1.2 `src/acmg_classifier/utils/progress.py` — モジュール全体 (H)
- `make_progress()` / `progress_iter()` がどこからも import されていない。
- 現在の pipeline は structlog でログを出しているのみで rich プログレスバーは未使用。
- 必要なら再追加すれば良いので、現状は削除候補。

### 1.3 `src/acmg_classifier/criteria/pathogenic/pp5.py` — `PP5Evaluator` クラス (M)
- `CriteriaRegistry._build_evaluators()` で意図的に登録されていない（ClinGen SVI が PP5/BP6 の廃止を勧告したため）。
- 自動評価ループでは絶対に呼ばれない。
- テスト (`tests/unit/test_criteria_pathogenic.py`) からのみ参照されている。
- 互換性のために残すなら現状維持。完全削除するなら enum (`ACMGCriterion.PP5`) と DEFAULT_STRENGTH も併せて要検討。
- マニュアルサプリメント経由での PP5 注入は `manual.py` 系では対応していない（manual list にも入っていない）。

### 1.4 `src/acmg_classifier/cli.py` — `_make_config()` (H)
- `cli.py` 内で定義のみ、どこからも呼ばれていない。
- 各サブコマンドが `Config(...)` を直接構築しているため死コード。

### 1.5 `src/acmg_classifier/classification/classifier_2015.py` — `_count()` / `_has()` (M)
- どちらもファイル内のみで定義されているが、`Classifier2015.classify()` から呼び出されていない。
- テストで使われている可能性があるので、削除前に `grep -r "_count\|_has"` で確認推奨。

---

## 2. 到達不能 (dead branch) / 効果のないコード

### 2.1 `src/acmg_classifier/pvs1/decision_tree.py:204` — `exon_skip_possible = True` ハードコード (H)
```python
exon_skip_possible = True  # conservative assumption without RNA data
if exon_skip_possible:
    ...
return (CriterionStrength.SUPPORTING, ...)  # ← 到達不可
```
- `exon_skip_possible` が常に `True` のため、最後の `return SUPPORTING` ブロックは到達しない。
- 2 つの選択肢:
  - (a) RNA データがあるときの分岐を実装する（変数を本物の判定に置き換える）。
  - (b) 変数と到達不能 return を削除して `Supporting` フォールバックは無いと明示する。
- スコアに影響するため、ドメイン判断必須。

---

## 3. 渡しても無視されるパラメータ

### 3.1 `src/acmg_classifier/local_db/vep_runner.py:137` — `amino_acid_change=aa_change` (H)
- `ConsequenceInfo` モデルには `amino_acid_change` フィールドが存在しない（`models/annotation.py:27` 参照）。
- Pydantic v2 のデフォルト `extra="ignore"` で黙って破棄される。
- すぐ上の `aa_change` 構築コード (3 行) も丸ごと死コード。
- 修正案:
  - (a) `ConsequenceInfo` に `amino_acid_change: Optional[str] = None` を追加して活用する。
  - (b) `aa_change` 構築と引数渡しの両方を削除する。

### 3.2 `src/acmg_classifier/models/criteria.py` — `not_met()` の `default_strength` 変数 (L)
```python
default_strength, direction = DEFAULT_STRENGTH[criterion]
return cls(..., strength=CriterionStrength.NOT_MET, direction=direction, ...)
```
- `default_strength` は読み出すが使わずに捨てている。
- `_, direction = DEFAULT_STRENGTH[criterion]` でよい（マイクロ最適化レベル）。

---

## 4. 重複ロジック / 統合候補

### 4.1 PP3 / BP4 の missense 分岐コード (M)
- `criteria/pathogenic/pp3.py` と `criteria/benign/bp4.py` の missense 分岐がほぼ対称形。
  - SpliceAI 抑制チェック → ESM1b / AlphaMissense 分岐の流れが全く同じ構造。
- 共通ヘルパー関数 `_score_missense(annotation, cfg, direction)` などに抽出すると重複が減る。
- 動作は変わらないので「リファクタリング」扱い。スコア計算の責任分界点が崩れる可能性もあるので慎重に。

### 4.2 ManualPathogenicEvaluator / ManualBenignEvaluator (L)
- `criteria/pathogenic/manual.py` と `criteria/benign/manual.py` のロジック本体（supplement から met / not_met を構築）は完全に同じ。
- `_MANUAL_CRITERIA` タプルだけが異なる。
- `BaseManualEvaluator(criteria_tuple)` のような共通基底に抽出可能。

### 4.3 `_load_thresholds` の TSV パース (L)
- `criteria/benign/bs1.py` 内の TSV パースが、`io/supplement_reader.py` 等で行っているパースと類似している。
- 直接の共通化メリットは小さい（カラム名・エラーポリシーが違う）。

---

## 5. 表記揺れ / 軽微なクリーンアップ

### 5.1 `criteria/registry.py` — PP5 と PP4 のコメント整合性 (L)
- 33 行目に `# PP5 (reputable-source) intentionally NOT registered` と書かれており説明済み。
- PP4 はマニュアル評価器 (`ManualPathogenicEvaluator`) に含まれているが、PP5 は完全に除外という非対称が分かりにくい。
- README へのリンクを足せば十分。

### 5.2 `models/classification.py:28` — `annotation: Optional[object]` 型注釈 (L)
- 「循環 import 回避のため」というコメント付きだが、Pydantic v2 では `from __future__ import annotations` と `TYPE_CHECKING` import で正しく型付け可能。
- 現状でも動作はする。

---

## 6. 参考: ドキュメント / コメント上の指摘

- `pvs1/nmd_predictor.py` の `predicts_nmd` — 「penultimate exon の最後 50bp」のチェックは未実装で、保守的に NMD ありとしている。コメントには書かれているが、将来 RNA データを使う場合の TODO 候補。
- `local_db/splice/squirls_predictor.py` — SQUIRLS 閾値は「approximate」とコメント済みだが、Walker calibration が来た時点で差し替えるべき箇所。

---

## 削除しなかった理由 (補足)

- `criteria/pathogenic/pp5.py` を残した理由: enum / strength マップ / 互換性のため。完全削除には enum 変更を伴う破壊的変更が必要。

