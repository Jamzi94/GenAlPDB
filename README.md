# GenAlPDB

> **AlphaFold 3 JSON & SLURM workflow** — generate AF3 input JSON from UniProt IDs, optionally fuse MSA, submit jobs to SLURM, parse runtimes, and plot performance.

<p align="center">
  <a href="#-overview--description-en">English</a> •
  <a href="#-aperçu--description-fr">Français</a>
</p>

---

## 🇬🇧 Overview / Description (EN)
GenAlPDB is an advanced pipeline adapted to an AlphaFold 3 workflow. For any set of UniProt IDs, it generates AF3 JSON, optionally fuses MSA data, submits jobs to SLURM, parses runtime metrics, and produces performance plots.

### Project Structure
```
dignose_workflow.py      : Diagnostic of the AF3 workflow (validates folders/files, script coherence, JSON health, recommendations).
remote_job_submit.py     : AF3 SLURM job submission (detects JSON, groups by base_id, sets MSA→inference dependencies).
slurm_job_submit.sh      : SLURM shell submit script (reads `config/benchmark_config.csv`, creates logs/output/msa_cache, launches MSA then inference).

scripts/submit_benchmark.py : High-level entry point: runs `scripts/main.py` (checks templates, (re)generates JSON, calls submission).
scripts/main.py             : Orchestrator (CSV template checks, JSON generation, SLURM hooks).
scripts/generate_af3_json.py: Builds AF3 input JSONs from `config/af3_template.csv` (fetches UniProt sequences, handles ligand/ion multiplicities, normalizes via `ligand_resolver`).
scripts/fusion_json.py      : Fuses an MSA JSON into the original AF3 JSON (validates A–Z chain IDs, deep recursive merge).
scripts/parse_slurm_time.py : Parses SLURM logs → `logs/runtime_summary.csv` (+ `.xlsx` if available).
scripts/plot_curves.py      : Performance plots (GPU ↔ time/memory) from the runtime summary (PNG/PDF).
scripts/ligand_resolver.py  : Ligand/Ion resolution: names → SMILES via `config/smiles_database.csv`, supports CCD codes.
scripts/generate_excel.py   : Tabular post-processing utilities (in progress).

config/af3_template.csv         : Input template (columns `uniprot_id`, `ligands`, `ions`).
config/af3_template_results.csv : Example results (`uniprot_id`, `gene`, `status`, `ligands`, `ions`, `json_file`).
config/benchmark_config.csv     : HPC resources (`gpu`, `cpus`, `mem`, `nodes`, `partition`, `walltime`).
config/smiles_database.csv      : Local SMILES database (`compound_name`, `smiles`, `source_sheet`).

README.md (root) : Original documentation, integrated into this adapted README.
```

### Installation
Clone the repository, (optionally) create a virtualenv, and install dependencies:
```bash
git clone <your-repo-url>
cd <your-repo>
python -m venv .venv && source .venv/bin/activate  # Windows: .venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### Pipeline Workflow
**Data Retrieval & JSON Generation (`scripts/generate_af3_json.py`)**
```
Reads `config/af3_template.csv` (columns: `uniprot_id`, `ligands`, `ions`).

Examples:
  ligands:  "ADP:2"           (single ligand with multiplicity)
            "ADP:2;ATP:2"     (multiple ligands separated by `;`)
  ions:     "MG:2"            (multiple ions supported with `;`, e.g. "MG:2;NA:1")

Fetches UniProt sequences, expands multiplicities, resolves ligands/ions via `ligand_resolver`
and `config/smiles_database.csv`. Writes AF3 JSON into `af3_json/` (e.g., `ABCC1_P33527_ATP2_MG2.json`).
```

**Optional MSA Fusion (`scripts/fusion_json.py`)**
```
Merges MSA-enriched JSON into the original AF3 JSON.
Validates/fixes chain IDs (A–Z), performs deep merge to preserve all fields.
```

**SLURM Submission (`remote_job_submit.py`, `slurm_job_submit.sh`)**
```
`remote_job_submit.py` scans `af3_json/`, groups JSON by base protein ID, and calls `slurm_job_submit.sh`.
The shell script reads `config/benchmark_config.csv`, creates `logs/`, `output/`, `msa_cache/`,
and submits MSA then inference jobs with proper dependencies.
```

**Runtime Parsing (`scripts/parse_slurm_time.py`)**
```
Parses SLURM logs and exports `logs/runtime_summary.csv` (+ `.xlsx` if supported).
```

**Visualization (`scripts/plot_curves.py`)**
```
Uses the runtime summary to generate performance plots (`time_vs_gpu.png/.pdf`, `memory_vs_gpu.png`).
```

**Diagnostics (`dignose_workflow.py`)**
```
Validates directory layout, presence of configs/templates, JSON health, and script consistency; prints actionable fixes.
```

### Usage
Prepare inputs in `config/`:
```
- af3_template.csv:      columns `uniprot_id`, `ligands`, `ions`
- benchmark_config.csv:  columns `gpu`, `cpus`, `mem`, `nodes`, `partition`, `walltime`
- smiles_database.csv:   columns `compound_name (=Simplified Molecular-Input Line-Entry System or SMILES)`, `smiles`, `source_sheet`
```

Run the full workflow (non-SLURM orchestrator):
```bash
python scripts/submit_benchmark.py
```

Run remote submission directly (optional):
```bash
python remote_job_submit.py [<project_root>]
bash slurm_job_submit.sh [<project_root>]
```

Generate JSON only:
```bash
python scripts/generate_af3_json.py
```

Fuse MSA JSON:
```bash
python scripts/fusion_json.py <json_original> <json_msa>
```

Parse runtimes:
```bash
python scripts/parse_slurm_time.py
```

Plot performance:
```bash
python scripts/plot_curves.py
```

### Options
```
remote_job_submit.py positional:  [root_dir]                 — project root (defaults to CWD)
fusion_json.py positional:        <json_original> <json_msa> — files to merge
-v:                               Verbose logging (if enabled in your logging setup)
```

### Output
```
af3_json/  : Generated AF3 input JSON files.
logs/      : Runtime logs and `runtime_summary.csv` (+ `.xlsx` when supported).
output/    : SLURM job outputs/artifacts.
msa_cache/ : Cache for MSA runs (when applicable).
```

> All sheets will include clickable hyperlinks (UniProt, RCSB, AlphaFold, PubMed…) where applicable.

### Notes
```
- UniProt API access may be required during JSON generation.
- SLURM wrapper is provided (`remote_job_submit.py` + `slurm_job_submit.sh`).
- Root-level files are fully integrated to reflect the actual repository layout.
```

---

## 🇫🇷 Aperçu / Description (FR)
GenAlPDB est un pipeline avancé adapté au workflow AlphaFold 3. À partir d’IDs UniProt, il génère des JSON AF3, fusionne optionnellement la MSA, soumet les jobs à SLURM, extrait les métriques d’exécution et produit des courbes de performance.

### Structure du projet
```
dignose_workflow.py      : Diagnostic complet du workflow AF3 (validation dossiers/fichiers, cohérence scripts, état des JSON, recommandations).
remote_job_submit.py     : Soumission des jobs AF3 à SLURM (détection JSON, groupement par base_id, dépendances MSA→inférence).
slurm_job_submit.sh      : Script shell de soumission SLURM (lit `config/benchmark_config.csv`, crée logs/output/msa_cache, lance MSA puis inférence).

scripts/submit_benchmark.py : Point d’entrée haut-niveau : lance `scripts/main.py` (vérifie templates, (re)génère JSON, appelle la soumission).
scripts/main.py             : Orchestration (vérification des CSV modèles, génération JSON, hooks SLURM).
scripts/generate_af3_json.py: Génère les JSON d’entrée AF3 depuis `config/af3_template.csv` (récupère séquences UniProt, gère multiplicité ligands/ions, normalise via `ligand_resolver`).
scripts/fusion_json.py      : Fusionne un JSON MSA dans le JSON AF3 original (validation ID de chaîne A–Z, fusion récursive robuste).
scripts/parse_slurm_time.py : Extrait les temps d’exécution à partir des logs SLURM → `logs/runtime_summary.csv` (+ `.xlsx` si possible).
scripts/plot_curves.py      : Courbes perfs (GPU ↔ temps/mémoire) à partir du résumé d’exécution (PNG/PDF).
scripts/ligand_resolver.py  : Résolution ligands/ions : noms → SMILES via `config/smiles_database.csv`, prise en charge des codes CCD.
scripts/generate_excel.py   : Utilitaires de post-traitement tabulaire (en cours).

config/af3_template.csv         : Modèle d’entrée (colonnes `uniprot_id`, `ligands`, `ions`).
config/af3_template_results.csv : Exemple de résultats (`uniprot_id`, `gene`, `status`, `ligands`, `ions`, `json_file`).
config/benchmark_config.csv     : Ressources HPC (`gpu`, `cpus`, `mem`, `nodes`, `partition`, `walltime`).
config/smiles_database.csv      : Base SMILES locale (`compound_name`, `smiles`, `source_sheet`).

README.md (racine) : Documentation fournie, intégrée dans ce README adapté.
```

### Installation
Cloner le dépôt, (optionnel) créer un environnement virtuel, puis installer les dépendances :
```bash
git clone <url-du-dépôt>
cd <votre-répertoire>
python -m venv .venv && source .venv/bin/activate  # Windows : .venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### Flux du pipeline
**Récupération & Génération JSON (`scripts/generate_af3_json.py`)**
```
Lit `config/af3_template.csv` (colonnes : `uniprot_id`, `ligands`, `ions`).

Exemples :
  ligands :  "ADP:2"
             "ADP:2;ATP:2"      (plusieurs ligands séparés par `;`)
  ions :    "MG:2"              (plusieurs ions possibles avec `;`, ex. "MG:2;NA:1")

Récupère les séquences UniProt, gère les multiplicités, normalise ligands/ions via `ligand_resolver`
et `config/smiles_database.csv`. Écrit des JSON AF3 dans `af3_json/` (ex. `ABCC1_P33527_ATP2_MG2.json`).
```

**Fusion MSA (optionnelle) (`scripts/fusion_json.py`)**
```
Fusionne le JSON enrichi MSA dans le JSON AF3 d’origine.
Valide/corrige les IDs de chaîne (A–Z), fusion récursive conservant tous les champs.
```

**Soumission SLURM (`remote_job_submit.py`, `slurm_job_submit.sh`)**
```
`remote_job_submit.py` parcourt `af3_json/`, groupe les JSON par ID protéique de base, et appelle `slurm_job_submit.sh`.
Le script shell lit `config/benchmark_config.csv`, crée `logs/`, `output/`, `msa_cache/`,
et soumet MSA puis inférence avec dépendances.
```

**Parsing des temps (`scripts/parse_slurm_time.py`)**
```
Analyse les logs SLURM et exporte `logs/runtime_summary.csv` (+ `.xlsx` si disponible).
```

**Visualisation (`scripts/plot_curves.py`)**
```
Génère des figures de performance (`time_vs_gpu.png/.pdf`, `memory_vs_gpu.png`).
```

**Diagnostic (`dignose_workflow.py`)**
```
Valide l’arborescence, la présence des configs/templates, l’état des JSON et la cohérence des scripts ; propose des correctifs actionnables.
```

### Utilisation
Préparer les entrées dans `config/` :
```
- af3_template.csv :      colonnes `uniprot_id`, `ligands`, `ions`
- benchmark_config.csv :  colonnes `gpu`, `cpus`, `mem`, `nodes`, `partition`, `walltime`
- smiles_database.csv :   colonnes `compound_name (=Simplified Molecular-Input Line-Entry System or SMILES)`, `smiles`, `source_sheet`
```

Lancer le workflow complet :
```bash
python scripts/submit_benchmark.py
```

Soumission distante directe (optionnel) :
```bash
python remote_job_submit.py [<racine_projet>]
bash slurm_job_submit.sh [<racine_projet>]
```

Générer uniquement les JSON :
```bash
python scripts/generate_af3_json.py
```

Fusionner un JSON MSA :
```bash
python scripts/fusion_json.py <json_original> <json_msa>
```

Parser les temps :
```bash
python scripts/parse_slurm_time.py
```

Tracer les perfs :
```bash
python scripts/plot_curves.py
```

### Options
```
remote_job_submit.py positionnel :  [root_dir]                 — racine du projet (par défaut : CWD)
fusion_json.py positionnels   :     <json_original> <json_msa> — fichiers à fusionner
-v :                                Logs verbeux (si activé dans votre configuration de logging)
```

### Sorties
```
af3_json/  : JSON d’entrée AF3 générés.
logs/      : Journaux d’exécution et `runtime_summary.csv` (+ `.xlsx` si disponible).
output/    : Sorties/artéfacts des jobs SLURM.
msa_cache/ : Cache pour les exécutions MSA (si pertinent).
```

> Les feuilles de résultats contiennent des hyperliens cliquables (UniProt, RCSB, AlphaFold, PubMed…) lorsque pertinent.

### Notes
```
- L’accès à l’API UniProt peut être requis lors de la génération des JSON.
- Le wrapper SLURM est fourni (`remote_job_submit.py` + `slurm_job_submit.sh`).
- Les fichiers situés à la racine sont pleinement intégrés pour refléter fidèlement la structure du dépôt.
```
