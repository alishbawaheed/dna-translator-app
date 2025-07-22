import streamlit as st
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight
from collections import Counter
import matplotlib.pyplot as plt
import re
import io

# ---------------- Utilities ---------------- #

def clean_dna(raw_text: str) -> str:
    """
    Take raw pasted FASTA or plain sequence text and return
    an uppercase string containing only A,T,G,C.
    """
    if not raw_text:
        return ""
    # Remove FASTA header lines starting with >
    lines = raw_text.strip().splitlines()
    if lines and lines[0].startswith(">"):
        lines = lines[1:]
    seq = "".join(lines)
    # Keep only ATGC (case-insensitive)
    seq = re.sub(r"[^ATGCatgc]", "", seq)
    return seq.upper()

def translate_dna(sequence: str, frame: int = 0) -> str:
    """
    Translate DNA starting at a given frame (0,1,2).
    Returns string of amino acids using Biopython table 1.
    """
    if frame not in (0,1,2):
        raise ValueError("frame must be 0, 1, or 2")
    trimmed = sequence[frame:]
    # Length must be multiple of 3 for clean translation; truncate remainder
    codon_len = len(trimmed) - (len(trimmed) % 3)
    coding = trimmed[:codon_len]
    if not coding:
        return ""
    # Use Biopython; it'll error on invalid chars, so we assume cleaned input
    return str(Seq(coding).translate(to_stop=False))

def find_start_stop_codons(sequence):
    starts, stops = [], []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if codon == "ATG":
            starts.append(i)
        elif codon in ["TAA", "TAG", "TGA"]:
            stops.append(i)
    return starts, stops

def find_orfs(sequence):
    """
    Return list of tuples: (start, end, length, frame)
    end is index AFTER stop codon (python slicing convention)
    """
    orfs = []
    for frame in range(3):
        in_orf = False
        start = None
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if not in_orf and codon == "ATG":
                in_orf = True
                start = i
            elif in_orf and codon in ["TAA", "TAG", "TGA"]:
                end = i + 3
                orfs.append((start, end, end - start, frame + 1))
                in_orf = False
                start = None
        # If sequence ends while in ORF and no stop codon, we ignore (partial ORF)
    return orfs

def compare_sequences(seq1, seq2):
    muts = []
    for i in range(min(len(seq1), len(seq2))):
        if seq1[i] != seq2[i]:
            muts.append((i+1, seq1[i], seq2[i]))
    return muts

def codon_usage(sequence):
    counts = Counter(sequence[i:i+3] for i in range(0, len(sequence)-2, 3))
    # Remove partial codon at end automatically since we stop at len-2
    return dict(counts)

def is_valid_protein(seq):
    valid = set("ACDEFGHIKLMNPQRSTVWY")
    return "".join([aa for aa in seq if aa in valid])

def guess_protein_name(protein_seq):
    patterns = {
        "MKT": "Kinase-like",
        "MEL": "Elongation Factor-like",
        "MAK": "Actin-like",
        "MGA": "G-protein-like",
        "MEN": "Enzyme-like",
        "M": "Generic Protein"
    }
    for pat, name in patterns.items():
        if protein_seq.startswith(pat):
            return name
    return "Unknown Protein"

def plot_codon_usage(usage_dict):
    if not usage_dict:
        st.info("No codons to display.")
        return
    codons = list(usage_dict.keys())
    counts = [usage_dict[c] for c in codons]
    fig, ax = plt.subplots()
    ax.bar(codons, counts)
    ax.set_xlabel("Codon")
    ax.set_ylabel("Count")
    ax.set_title("Codon Usage")
    plt.xticks(rotation=90)
    st.pyplot(fig)

def fasta_from_protein(protein_seq, name="Translated_Protein"):
    wrapped = "\n".join(protein_seq[i:i+60] for i in range(0, len(protein_seq), 60))
    return f">{name}\n{wrapped}\n"

# ---------------- Streamlit App ---------------- #

st.title("🧬 DNA Translator & ORF Finder")

raw_seq = st.text_area("Paste DNA or FASTA sequence:", height=200)
dna_seq = clean_dna(raw_seq)

if raw_seq and not dna_seq:
    st.warning("No valid A/T/G/C bases found after cleaning input. Please check your sequence.")

if dna_seq:
    # RNA
    st.subheader("🧬 RNA Sequence")
    rna_seq = dna_seq.replace("T", "U")
    st.code(rna_seq)

    # Start/Stop Codons
    st.subheader("🧬 Start & Stop Codons")
    start_pos, stop_pos = find_start_stop_codons(dna_seq)
    st.write(f"Start Codons (ATG): {len(start_pos)} at positions: {start_pos}")
    st.write(f"Stop Codons (TAA/TAG/TGA): {len(stop_pos)} at positions: {stop_pos}")

    # Codons Used (list)
    st.subheader("🧬 Codons Used (list)")
    codon_list = [dna_seq[i:i+3] for i in range(0, len(dna_seq)-2, 3)]
    st.write(codon_list)

    # Codon Usage Chart
    st.subheader("📊 Codon Usage Chart")
    usage = codon_usage(dna_seq)
    plot_codon_usage(usage)

    # ORFs
    st.subheader("🧬 Open Reading Frames (ORFs)")
    orfs = find_orfs(dna_seq)
    if orfs:
        for idx, (start, end, length, frame) in enumerate(orfs, 1):
            st.write(f"ORF {idx}: Frame {frame}, Start={start}, Stop={end}, Length={length} bp")
    else:
        st.info("No complete ORFs found (ATG ... Stop).")

    if orfs:
        longest = max(orfs, key=lambda x: x[2])
        st.success(f"Longest ORF: Frame {longest[3]}, Start {longest[0]}, Stop {longest[1]}, Length {longest[2]} bp")

    # Protein Translations
    st.subheader("🧬 Protein Translations (3 Reading Frames)")
    for frame in range(3):
        raw_protein = translate_dna(dna_seq, frame)
        clean_protein = is_valid_protein(raw_protein)
        st.markdown(f"**Frame {frame+1}:**")
        st.code(clean_protein if clean_protein else "(no valid amino acids)")

        # Protein name guess
        st.markdown(f"**Protein Name Guess (Frame {frame+1}):** {guess_protein_name(clean_protein)}")

        # Molecular weight (only if we have valid AA)
        if clean_protein:
            mw = molecular_weight(clean_protein, seq_type='protein')
            st.markdown(f"**Molecular Weight:** {mw:.2f} Da")
        else:
            st.markdown("**Molecular Weight:** n/a (ambiguous or empty sequence)")

    # Mutation Detection
    st.subheader("🧬 Mutation Detector (Optional)")
    ref_raw = st.text_input("Paste another DNA sequence to compare (same length recommended):")
    if ref_raw:
        ref_seq = clean_dna(ref_raw)
        if not ref_seq:
            st.warning("No valid A/T/G/C bases in comparison sequence.")
        else:
            muts = compare_sequences(ref_seq, dna_seq)
            if muts:
                st.write(f"Found {len(muts)} mutation(s):")
                for pos, a, b in muts:
                    st.write(f"Position {pos}: {a} → {b}")
            else:
                st.success("No differences detected.")

    # FASTA Download (Frame 1 protein)
    st.subheader("⬇️ Download Protein Sequence (FASTA)")
    frame1_protein = is_valid_protein(translate_dna(dna_seq, 0))
    fasta_data = fasta_from_protein(frame1_protein, "Translated_Protein_Frame1")
    st.download_button("Download FASTA", data=fasta_data, file_name="protein_frame1.fasta", mime="text/plain")

    # Reverse Complement
    st.subheader("🔁 Reverse Complement")
    if st.button("Show Reverse Complement"):
        st.code(str(Seq(dna_seq).reverse_complement()))
else:
    st.info("Please paste a DNA or FASTA sequence above to begin.")

