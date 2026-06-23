# 残項目の実装要件（データソース / クエリ）

凡例: **[S]** 追加データ不要・実装可能 / **[D]** 外部データ要（ソース提示）/ **[Q]** 導出スクリプト要（コード提示・そちらで実行）

---

## A3. PM4 region / strength

ほぼ cspec JSON 内に座標・残基・強度があり、PhyloP は既存パイプライン（`PhyloPReader`）、欠失内容は ClinVar SQLite ランタイムクエリ（私が新規関数を実装）で対応可能。**大半 [S]**。

新スキーマ案（`disease_prevalence.tsv`）: `pm4_strong_residues`, `pm4_strong_regions`, `pm4_region`(allow-list), `pm4_deny_regions`, `pm4_stoploss`(strength), `pm4_conserved_phylop`(cutoff), `pm4_indel_only`/`pm4_stoploss_only`(applicability限定)。

| 遺伝子 | 内容 | 区分 | ソース/備考 |
|---|---|---|---|
| RUNX1 (GN008) | PM4_Strong/Mod 残基 R107,K110,A134,R162,R166,S167,R169,G170,K194,T196,D198,R201,R204；Supporting RHD 89-204；stop-loss→PVS1拡張 | [S] | 残基は JSON 内 |
| MYOC (GN019) | olfactomedin AA246-502；>10%→Mod / ≤10%→Sup（%は aa_len から算出） | [S] | JSON 内 |
| DICER1 (GN024) | RNase IIIb p.Y1682-S1846→Mod；外側→Sup；repeat除外 p.D606-609/E1418-1420/E1422-1425 | [S] | JSON 内 |
| RPGR (GN106) | stop-loss aa1153→Strong；exons1-14 + ORF15 aa585-1078→Mod | [S] | exon は MANE GFF（保有）から生成 |
| PTEN (GN003) | catalytic motif（PM1参照）+ protein extension | [S] | PM1 hotspot 流用 |
| CTLA4/PIK3R1/RPE65/GUCY2D/AIPL1/ABCA4 | ≥2aa かつ PhyloP保存 (CTLA4/PIK3R1≥2, ABCA4≥7.367)；SpliceAI<0.2 | [S] | 既存 phyloP bigwig + 既存 splice |
| CDKL5/FOXG1/MECP2 | deny-region（Pro-rich 等） | [S] | JSON 内（座標あり） |
| CDH1/ATM | stop-loss のみ適用 | [S] | JSON 内 |
| CYP1B1 | stop-loss を PM4 から除外 | [S] | JSON 内 |
| OTC/LDLR/RS1 | サイズ上限（<1 exon）/ PM2併用 caveat | [S] | exon は MANE |
| **SCID群** (FOXN1/ADA/DCLRE1C/IL7R/JAK3/RAG1/RAG2/IL2RG) | 欠失領域が既知P/LP含む→Mod、VUS含む→Sup | **[S]** | ClinVar SQLite に**新規ランタイムクエリ**を実装（`query_clinvar_in_protein_range(gene, lo, hi, sig)`）。外部データ不要 |
| VHL (GN078) | B/α/β2 ドメイン in-frame indel | **[D]** | JSON は「B and alpha domains」のみで**数値座標なし**。座標 β63-155/α156-192/β2 193-204 の出典 = VHL VCEP 公開資料（PMID:20560986 / cspec GN078 添付）。座標 3 行いただければ実装 |

---

## B1. GALT PVS1（in-frame exon 6/7/9 splice 例外）

- exon 座標表: **[S]** — `build_vcep_pvs1_exons.py` の `TARGET_TRANSCRIPTS` に `"GALT": "NM_000155"` を追加し、保有の MANE GFF から再生成可能。
- per-exon の強度（exon 6/7/9 skip → Strong/Moderate のどちら）: **[D]**
  - **ソース**: cspec GN158 の「PVS1 decision tree **attached**」PDF（JSON 本文に強度が無く "Refer to the attachment" のみ）。
  - **必要フォーマット**（`vcep_pvs1_splice_exons.tsv` に追記）:
    ```
    GALT	6	<very_strong|strong|moderate|supporting|na>	<note>
    GALT	7	...
    GALT	9	...
    ```
  - 添付の決定木 PDF（または exon6/7/9 の強度 3 値）をいただければ実装。

---

## B2. PM1（LDLR / VHL / FBN1）

| 遺伝子 | 必要データ | 区分 | ソース / 導出 |
|---|---|---|---|
| LDLR (GN013) | exon4 の codon 範囲 + **60 conserved Cys 残基** | exon4=[S]（MANE）／Cys=[D] or [Q] | Cys リストの正典 = LDLR VCEP 論文 **Supp. Table 4**。無ければ [Q]: UniProt **P01130** の LDL-receptor class A repeat 内 Cys を抽出（下記スクリプト）。ただし「highly conserved」サブセットは Supp Table 4 が正。|
| VHL (GN078) | Germline/Somatic hotspot 表 + cancerhotspots.org ≥10 instances | **[D]** | (1) VHL VCEP「Germline and Somatic Hotspots」表（cspec GN078 添付）。(2) cancerhotspots.org のダウンロード（https://www.cancerhotspots.org/#/download）。VHL の残基別 instance 数を CSV でいただければ ≥10/<10 を Mod/Sup に振り分け実装。|
| FBN1 (GN022) | cbEGF/EGF/TB/hybrid ドメイン境界 + ドメイン内 Cys 位置、Ca結合/水酸化残基、Gly motif | **[D] or [Q]** | UniProt **P35555** の domain/feature。motif ルール（Cys 種別・(D/N)-X-(D/N)... コンセンサス）は現エンジンの「残基/範囲」表現を超えるため、**ドメイン別 Cys 残基リスト**へ落としてもらえれば PM1 残基として実装可。下記スクリプトで UniProt から cbEGF ドメイン Cys を抽出可能。|

---

## C1. PS1 paralog / analogous-residue

| グループ | 区分 | 内容 |
|---|---|---|
| RASopathy 同番号（Group1 HRAS/NRAS/KRAS, Group2 MAP2K1/MAP2K2, Group3 SOS1/SOS2）+ HBA2/HBA1 | **[S]** | グループは JSON 内。これらは番号が一致するため、PS1 評価器を「同一 hgvs_p を兄弟遺伝子でも ClinVar 照会」に拡張すれば実装可。新規 `ps1_paralog_group` 列（例 `HRAS,NRAS,KRAS`）+ クエリ拡張。外部データ不要 |
| SCN1A/2A/3A/8A, KCNQ1 群 | **[Q]** | 番号が一致しないため **paralog 残基対応表**が必要。導出 = UniProt 配列の多重整列。下記スクリプトで対応表 CSV を生成（そちらで実行）→ `gene_a, pos_a, gene_b, pos_b` 形式でいただければ実装 |

---

## C2. 非REVEL in-silico cutoff（実装には新スコア源が必要）

| 予測器 | 対象 | 区分 | ソース |
|---|---|---|---|
| BayesDel | TP53(GN009), BRCA1(GN092), BRCA2(GN097) の PP3/BP4 | **[D]** | BayesDel 事前計算スコア（https://fengbj-laboratory.org/BayesDel/ ; addAF / noAF 版）。VCF/TSV で変異→スコアの参照テーブルが必要。パイプライン（`AnnotationData`）に `bayesdel` フィールド追加が前提 |
| HCI-prior (MAPP+PP2) | MLH1/MSH2/MSH6(GN115/137/138), PMS2(GN139) | **[N/A]** | InSiGHT/HCI の事前確率（gene-specific prior probability of pathogenicity）。`hci-priors.hci.utah.edu/PRIORS`（web専用・APIなし）および `hci-lovd.hci.utah.edu` の `*_priors` LOVD DB。**LOVD 規約でバルクダウンロード不可**（「do not download the LOVD package; obtain from the respective databases」）→ **自動アノテーション不可・実装見送り**。per-variant の web/DB 参照値を手動 supplement TSV（`SupplementEntry`: PP3, strength）で投入するのが唯一の対応経路 |
| CADD | CTLA4/PIK3CD/PIK3R1/BMPR2/ABCA4 の補助条件（REVELと2of3等） | **[D]** | CADD スコア（https://cadd.gs.washington.edu/ ; GRCh38 whole-genome SNV TSV.gz）。`cadd` フィールド追加が前提 |
| AlphaMissense | BMPR2 等の 2of3 条件 | **[S]** | 既にパイプラインに存在（`annotation.alphamissense`）。条件ロジックのみ実装可 |

各スコアの参照テーブル（chrom,pos,ref,alt → score）をいただければ、annotation orchestrator にフィールドを足して評価器を実装します。

---

## 導出スクリプト（[Q] 用・そちらで実行して結果 CSV を返却）

### (a) LDLR / FBN1 の ドメイン内 Cys 残基（UniProt feature から）
```python
# 要: requests。UniProt の feature(domain/disulfide) から指定ドメイン内 Cys 位置を抽出。
import requests
def domain_cys(acc, domain_keywords):
    j = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}.json").json()
    seq = j["sequence"]["value"]
    cys = {i+1 for i,a in enumerate(seq) if a == "C"}
    doms = [f for f in j.get("features",[])
            if f["type"] in ("Domain","Repeat")
            and any(k.lower() in (f.get("description","")).lower() for k in domain_keywords)]
    out=set()
    for f in doms:
        lo=int(f["location"]["start"]["value"]); hi=int(f["location"]["end"]["value"])
        out |= {c for c in cys if lo<=c<=hi}
    return sorted(out)
print("LDLR class-A Cys:", domain_cys("P01130", ["LDL-receptor class A"]))
print("FBN1 cbEGF Cys:", domain_cys("P35555", ["EGF-like; calcium-binding","cbEGF"]))
```
→ 返却: 遺伝子ごとの Cys 残基リスト（カンマ区切り int）。

### (b) SCN / KCNQ paralog 残基対応表（多重整列）
```python
# 要: biopython + 配列。MAFFT 等で整列し、整列カラム経由で残基番号を対応付け。
# 入力: 各遺伝子の UniProt 参照配列（SCN1A P35498 / SCN2A Q99250 / SCN3A Q9NY46 / SCN8A Q9UQD0）
# 出力 CSV: gene_a,pos_a,gene_b,pos_b （同一整列カラムに乗る残基ペア）
# （MAFFT 実行環境が無ければ配列 FASTA をいただければこちらで整列します）
```
→ 返却: `gene_a,pos_a,gene_b,pos_b` の対応表 CSV。

### (c) VHL cancerhotspots ≥10 instances
```text
cancerhotspots.org → Download → "Cancer Hotspots v2" (single residue) を取得し、
VHL 行のみ抽出して  residue, instance_count  を CSV で返却。
```
→ 返却: `residue,count`（≥10 を Moderate、<10 を Supporting に振り分け）。
