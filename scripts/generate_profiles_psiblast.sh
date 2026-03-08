#!/bin/bash
# generate PSSM profiles for control proteins using PSI-BLAST
# runs against the partial uniref90 database (first volumes)
# uses 3 iterations with e-value threshold 0.001

PSIBLAST="/storage/kiran-stuff/protein-folds/tools/ncbi-blast-2.16.0+/bin/psiblast"
BLASTDB="/storage/kiran-stuff/NABind/databases/uniref90_blast/uniref90"
CONTROLS_DIR="/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination/structures/pnabind_controls"
WORK_DIR="/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination/structures/graphrbf_controls"

PROTEINS="P0A9X9 P62244 P03023 P04637 P67809 P69441 U1A_crystal MS2_crystal LambdaRep_crystal"

for PROT in $PROTEINS; do
    echo "=== $PROT ==="
    PDIR="$WORK_DIR/$PROT"
    FASTA="$PDIR/${PROT}.fasta"
    PSSM_OUT="$PDIR/${PROT}_A.pssm"

    # check if fasta exists, if not create from sequence
    if [ ! -f "$FASTA" ]; then
        # extract sequence from PDB
        python3 -c "
import sys
res_dict = {'GLY': 'G', 'ALA': 'A', 'VAL': 'V', 'ILE': 'I', 'LEU': 'L',
            'PHE': 'F', 'PRO': 'P', 'MET': 'M', 'TRP': 'W', 'CYS': 'C',
            'SER': 'S', 'THR': 'T', 'ASN': 'N', 'GLN': 'Q', 'TYR': 'Y',
            'HIS': 'H', 'ASP': 'D', 'GLU': 'E', 'LYS': 'K', 'ARG': 'R'}
seen = set()
seq = []
first_chain = None
with open('$CONTROLS_DIR/$PROT.pdb') as f:
    for line in f:
        if line.startswith('ATOM'):
            rn = line[17:20].strip()
            ch = line[21]
            if first_chain is None: first_chain = ch
            if ch != first_chain: continue
            key = (ch, line[22:27].strip())
            if key not in seen and rn in res_dict:
                seen.add(key)
                seq.append(res_dict[rn])
print(f'>${PROT}')
print(''.join(seq))
" > "$FASTA"
    fi

    # skip if PSSM already exists and is non-trivial
    if [ -f "$PSSM_OUT" ] && [ $(wc -c < "$PSSM_OUT") -gt 5000 ]; then
        # check if it's a real PSSM (not our dummy one)
        if grep -q "0.00 0.00" "$PSSM_OUT" | head -1 && ! grep -q "position-specific" "$PSSM_OUT" 2>/dev/null; then
            echo "  using existing PSSM"
            continue
        fi
    fi

    echo "  running PSI-BLAST (3 iterations, 16 threads)..."
    $PSIBLAST -query "$FASTA" \
        -db "$BLASTDB" \
        -num_iterations 3 \
        -evalue 0.001 \
        -num_threads 16 \
        -out_ascii_pssm "$PSSM_OUT" \
        -save_pssm_after_last_round \
        -out /dev/null 2>&1 | tail -3

    if [ -f "$PSSM_OUT" ]; then
        echo "  PSSM generated: $(wc -l < "$PSSM_OUT") lines"
    else
        echo "  ERROR: PSSM generation failed"
    fi
done

echo "Done!"
