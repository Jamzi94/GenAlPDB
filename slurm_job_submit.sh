#!/usr/bin/env bash
# Exit immediately on errors and treat unset variables as errors
set -euo pipefail

# Define paths and directories
ROOT="${1:-$PWD}"
CONFIG="$ROOT/config/benchmark_config.csv"
JSON_DIR="$ROOT/af3_json"
LOG_DIR="$ROOT/logs"
OUT_DIR="$ROOT/output"
CACHE_DIR="$ROOT/msa_cache"

# Create necessary directories
mkdir -p "$LOG_DIR" "$OUT_DIR" "$CACHE_DIR"
# Read configuration lines (skip header)
mapfile -t RES < <(tail -n +2 "$CONFIG")

# Function to recursively find *_data.json files in the cache directory
find_data_json() {
  local cache_dir="$1"
  local base_id="$2"
  local data_json_files=()
  
  # Recursive search for all *_data.json files in the cache directory and its subdirectories
  while IFS= read -r file; do
    data_json_files+=("$file")
  done < <(find "$cache_dir/$base_id" -type f -name "*_data.json" 2>/dev/null)
  
  # If files were found, return the first one
  if [ ${#data_json_files[@]} -gt 0 ]; then
    echo "${data_json_files[0]}"
    return 0
  fi
  
  # No file found
  return 1
}

# Extract base_id from filename (first two parts separated by _)
extract_base_id() {
  local filename="$1"
  local base_name=$(basename "$filename" .json)
  echo "$base_name" | cut -d_ -f1-2
}

echo "Starting job submission process..."
echo "Looking for valid JSON input files..."

# First pass: collect all base_ids and their corresponding JSON files
declare -A BASE_ID_FILES
for json in "$JSON_DIR"/*.json; do
  [[ -f "$json" ]] || continue
  
  # Skip chain mapping files if they still exist
  if [[ "$json" == *"_chain_mapping.json" ]]; then
    echo "→ Skipping mapping file: $json (not a valid input)"
    continue
  fi
  
  # Extract base_id
  base_id=$(extract_base_id "$json")
  
  # Add to array of files for this base_id
  BASE_ID_FILES[$base_id]+=" $json"
done

echo "Found ${#BASE_ID_FILES[@]} unique protein base IDs to process"

# Second pass: process each base_id once for MSA
for base_id in "${!BASE_ID_FILES[@]}"; do
  echo "→ Processing base_id: $base_id"
  
  # Check if MSA already exists
  data_json=""
  if data_json=$(find_data_json "$CACHE_DIR" "$base_id"); then
    echo "→ MSA file found for base_id $base_id: $data_json"
    
    # Create a script to handle MSA fusion for all JSON files of this base_id
    fusion_script="$LOG_DIR/${base_id}_fusion.sh"
    cat > "$fusion_script" <<EOF
#!/bin/bash
# This script fuses MSA data into all JSON files for base_id $base_id

MSA_FILE="$data_json"
echo "Using MSA file: \$MSA_FILE"

# Process each JSON file for this base_id
for json_file in ${BASE_ID_FILES[$base_id]}; do
  echo "Fusing MSA data into \$json_file"
  python3 "$ROOT/scripts/fusion_json.py" "\$json_file" "\$MSA_FILE"
  
  # Verify fusion was successful by checking for unpairedMsa field
  if grep -q "unpairedMsa" "\$json_file"; then
    echo "Fusion successful: unpairedMsa field found in \$json_file"
  else
    echo "WARNING: Fusion may have failed - unpairedMsa field not found in \$json_file"
  fi
done
EOF
    chmod +x "$fusion_script"
    
    # Execute fusion script immediately for existing MSA
    echo "→ Executing fusion script for existing MSA"
    bash "$fusion_script"
    
    # No MSA job needed
    msa_job_id="0"
  else
    # Submit MSA job only once per base_id
    echo "→ No MSA found for base_id $base_id, submitting new MSA job"
    
    # Use the first JSON file for this base_id for MSA calculation
    first_json=$(echo ${BASE_ID_FILES[$base_id]} | tr ' ' '\n' | sort | head -1)
    
    # Create a script to handle MSA fusion after MSA generation
    fusion_script="$LOG_DIR/${base_id}_fusion.sh"
    cat > "$fusion_script" <<EOF
#!/bin/bash
# This script fuses MSA data into all JSON files for base_id $base_id

# Find the MSA file
MSA_FILE=\$(find "$CACHE_DIR/$base_id" -type f -name "*_data.json" | head -1)

if [ -z "\$MSA_FILE" ]; then
  echo "ERROR: No MSA file found for $base_id"
  exit 1
fi

echo "Using MSA file: \$MSA_FILE"

# Process each JSON file for this base_id
for json_file in ${BASE_ID_FILES[$base_id]}; do
  echo "Fusing MSA data into \$json_file"
  python3 "$ROOT/scripts/fusion_json.py" "\$json_file" "\$MSA_FILE"
  
  # Verify fusion was successful by checking for unpairedMsa field
  if grep -q "unpairedMsa" "\$json_file"; then
    echo "Fusion successful: unpairedMsa field found in \$json_file"
  else
    echo "WARNING: Fusion may have failed - unpairedMsa field not found in \$json_file"
  fi
done
EOF
    chmod +x "$fusion_script"
    
    # Create MSA job script
    msa_slurm="$LOG_DIR/${base_id}_msa.slurm"
    
    # Get configuration from first entry in RES array
    IFS=',' read -r gpu cpus mem nodes partition walltime <<<"${RES[0]}"
    
    cat > "$msa_slurm" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${base_id}_msa
#SBATCH --partition=cpu,gpu-h100,gpu-l40,gpu-a40 
#SBATCH --gres=gpu:0
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mem=${mem}
#SBATCH --nodes=1
#SBATCH --time=${walltime}
#SBATCH --output=${LOG_DIR}/${base_id}_msa.out
#SBATCH --error=${LOG_DIR}/${base_id}_msa.err

# Load Conda
eval "$(/nfs/home/gomes26/miniconda3/bin/conda shell.bash hook)"
conda activate af3_benchmark

# Print debug info
echo "===== JOB STARTED ====="
echo "Running on node: $(hostname)"
echo "GPU status:"
nvidia-smi
echo "===== STARTING AF3 RUN ====="
module load alphafold/3
export JAX_PLATFORMS=cpu

# Run AlphaFold MSA generation
apptainer exec --nv \\
  \${AF_IMAGE_DIR}/alphafold3.sif \\
  python \${AF_BIN_DIR}/run_alphafold.py \\
    --model_dir=\${AF_MODEL_PARAMETERS_DIR} \\
    --db_dir=\${AF_DATABASES_DIR} \\
    --json_path=${first_json} \\
    --output_dir=${CACHE_DIR}/${base_id} \\
    --run_inference=false

# After MSA generation, run the fusion script
bash ${fusion_script}
EOF
    sed -i 's/\r$//' "$msa_slurm" && chmod +x "$msa_slurm"
    msa_job_id=$(sbatch --parsable "$msa_slurm")
    echo "    MSA job submitted (JobID=${msa_job_id})"
  fi
  
  # Third pass: process each JSON file for inference
  for json in ${BASE_ID_FILES[$base_id]}; do
    full_id=$(basename "$json" .json)
    
    echo "→ Processing inference for file: $json (base_id: $base_id)"
    
    for r in "${RES[@]}"; do
      IFS=',' read -r gpu cpus mem nodes partition walltime <<<"$r"
      job="${full_id}_${gpu}_${cpus}_${mem}_${nodes}_${partition}"
      infer_slurm="$LOG_DIR/${job}_infer.slurm"
  
      # Define dependency only if an MSA job was submitted and it's not a dummy ID
      dependency_line=""
      if [[ -n "${msa_job_id:-}" && "${msa_job_id}" != "0" ]]; then
        dependency_line="#SBATCH --dependency=afterok:${msa_job_id}"
      fi
      
      cat > "$infer_slurm" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${full_id}_infer
#SBATCH --partition=gpu-h100,gpu-l40,gpu-a40
#SBATCH --gres=gpu:${gpu}
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mem=20GB
#SBATCH --nodes=${nodes}
#SBATCH --time=${walltime}
${dependency_line}
#SBATCH --output=${LOG_DIR}/${full_id}.out
#SBATCH --error=${LOG_DIR}/${full_id}.err


# Load Conda
eval "$(/nfs/home/gomes26/miniconda3/bin/conda shell.bash hook)"
conda activate af3_benchmark

# Print debug info
echo "===== JOB STARTED ====="
echo "Running on node: $(hostname)"
echo "GPU status:"
nvidia-smi
echo "===== STARTING AF3 RUN ====="
module load alphafold/3

# Verify that the JSON file contains the unpairedMsa field
if ! grep -q "unpairedMsa" "${json}"; then
  echo "ERROR: The JSON file does not contain unpairedMsa field. MSA fusion may have failed."
  echo "Attempting emergency fusion..."
  
  # Find MSA file
  MSA_FILE=\$(find "${CACHE_DIR}/${base_id}" -type f -name "*_data.json" | head -1)
  
  if [ -n "\$MSA_FILE" ]; then
    echo "Found MSA file: \$MSA_FILE"
    echo "Running emergency fusion..."
    python3 "${ROOT}/scripts/fusion_json.py" "${json}" "\$MSA_FILE"
    
    if ! grep -q "unpairedMsa" "${json}"; then
      echo "Emergency fusion failed. Cannot proceed with inference."
      exit 1
    else
      echo "Emergency fusion successful. Proceeding with inference."
    fi
  else
    echo "No MSA file found. Cannot proceed with inference."
    exit 1
  fi
fi

apptainer exec --nv \\
  \${AF_IMAGE_DIR}/alphafold3.sif \\
  python \${AF_BIN_DIR}/run_alphafold.py \\
    --model_dir=\${AF_MODEL_PARAMETERS_DIR} \\
    --db_dir=\${AF_DATABASES_DIR} \\
    --json_path=${json} \\
    --output_dir=$OUT_DIR/${job}_out \\
    --run_data_pipeline=false
EOF
      sed -i 's/\r$//' "$infer_slurm" && chmod +x "$infer_slurm"
      sbatch "$infer_slurm"
      echo "    Inference job submitted (depends on MSA: ${msa_job_id:-none})"
    done
  done
done

echo "All jobs submitted."

