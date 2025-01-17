import streamlit as st
# --- Page Configurations ---
st.set_page_config(page_title="MultiClassClassificationInput App", layout="wide")

from backend import user_input_menu, multi_user_input_menu

# --- Sidebar Navigation ---
st.sidebar.title("Navigation")
pages = ["Home", "Start Task","Meta Data", "Glossary","Demonstration", "About"]
selected_page = st.sidebar.selectbox("Select Page :", pages)

# --- Home Page ---
if selected_page == "Home":
    st.title("Welcome to the Multi Class Classification App")
    st.write("**Home Page**")

    st.write("...")
    st.write("...")
    st.write("...")

# --- Start Task Page ---
elif selected_page == "Start Task":
    st.title("Start Task")
    st.write("**Begin the task by interacting with the backend process.**")
    col1,col2 = st.columns(2)

    with col1:
        tid = st.text_input("Enter the Transcript ID: ", placeholder="e.g., Ca_00001", key="Tid_input1").strip()

    with col2:
        mtid = st.text_input("Enter multiple Transcript IDs: ", placeholder="e.g., Ca_00001, Ca_00002", key="mTid_input2").strip()
        if mtid:
            mtid_list = [item.strip() for item in mtid.replace(",", " ").split()]
            mtid_list = list(set(mtid_list))
            mtid = ",".join(mtid_list)

    if st.button("Start"):
        if tid:
            result = user_input_menu(tid)
            st.write(result)
            st.toast("Task completed successfully.")
        elif mtid:
            result =multi_user_input_menu(mtid)
            st.write(result)
            st.toast("Task completed successfully.")
        else:
            st.warning("Need either a Transcript ID or multiple IDs to proceed.")
    elif tid == "":
        st.warning("Need Transcript ID to proceed.")
    else:
        st.write("Press the 'Start' button to begin the task.")
        st.write("Follow the instructions or check out demonstrations")

# --- Meta Data Page ---
elif selected_page == "Meta Data":
    st.title("Meta Data")
    st.write("**Key Insights and Analytics from the Application Backend**")
    col1, col2 = st.columns([1, 2])

    with col1:
        st.image("Images/1.png", caption="Expression Data Heatmap", use_container_width=True)
        st.write("")
        #st.markdown("### Another Image Title")
        st.write("")
        st.image("Images/2.png", caption="SVM Kernel performance", use_container_width=True)
        st.write("")
        st.write("")
        st.image("Images/7.png", caption="Tissue Specific Distribution plots", use_container_width=True)
        st.write("")
        st.write("")
        st.image("Images/5.png", caption="WGCNA heatmaps", use_container_width=True)

    with col2:
        st.image("Images/4.png", caption="Functional Annotation [Root Tissues]", use_container_width=True)
        st.write("")
        st.image("Images/6.png", caption="Comaprison of lncRNAs, TF and Non-TF", use_container_width=True)
    st.image("Images/3.png", caption="Performance charts for all files", use_container_width=True)

# --- Glossary Page ---
elif selected_page == "Glossary":
    st.title("Glossary")
    st.write("**Key Terms and Definitions**")
    glossary_entries = {
    'GO - Gene Ontology': '- a framework for the model of biology that describes gene functions in a species-independent manner.',
    'KEGG - Kyoto Encyclopedia of Genes and Genomes': '- a database resource for understanding high-level functions and utilities of biological systems.',
    'FPKM - Fragments Per Kilobase of transcript per Million mapped reads': '- a normalized method for counting RNA-seq reads.',
    'miRNA - MicroRNA': '- small non-coding RNA molecules that regulate gene expression by binding to complementary sequences on target mRNA.',
    'lncRNA - Long Non-Coding RNA': '- a type of RNA molecule that is greater than 200 nucleotides in length but does not encode proteins.',
    'ST - Seed Tissue': '- the tissue in seeds that supports the development of the embryo and storage of nutrients.',
    'FDS - Flower Development Stages': '- the various phases of growth and development that a flower undergoes from bud to bloom.',
    'FP - Flower Parts': '- the various components that make up a flower, including petals, sepals, stamens, and carpels.',
    'GT - Green Tissues': ' - plant tissues that are photosynthetic, primarily found in leaves and stems.',
    'RT - Root Tissues': '- the tissues found in the root system of a plant, involved in nutrient absorption and anchorage.',
    'TF - Transcription Factor': '- a protein that controls the rate of transcription of genetic information from DNA to messenger RNA.',
    'Non-TF - Non-Transcription Factors': '- proteins or molecules that do not directly bind to DNA to initiate or regulate transcription, but still influence gene expression through other mechanisms.',
    'WGCNA - Weighted Gene Co-expression Network Analysis': '- a method for finding clusters (modules) of highly correlated genes and studying their relationships to clinical traits.',
    'PPI - Protein-Protein Interaction': '- physical contacts between two or more proteins that occur in a living organism and are essential for various biological functions, including signal transduction and gene regulation.',
    'SNP CALLING - Single Nucleotide Polymorphism': 'The process of identifying single nucleotide polymorphisms (SNPs) in a genome from sequencing data. SNPs are variations at a single position in the DNA sequence, and SNP calling is crucial for genetic studies and disease association analyses.',
    'PEPTIDE SEQUENCE': 'A sequence of amino acids that make up a peptide, which is a short chain of amino acids linked by peptide bonds.',
    'CDS SEQUENCE - Coding Sequence': '- the portion of a gene\'s DNA or RNA that codes for a protein.',
    'TRANSCRIPT SEQUENCE': 'The RNA sequence transcribed from a gene, which may be translated into a protein or may function as non-coding RNA.',
    'GENOMIC SEQUENCE': 'The complete sequence of nucleotides (DNA or RNA) that make up the entire genome of an organism.'}

    for term, definition in glossary_entries.items():
        with st.expander(term):
            st.write(definition)
    
# --- Demonstration Page ---
elif selected_page == "Demonstration":
    st.title("Demonstration Page")
    st.write("**Learn how to use this interface**")

    # Add help content here
    st.write("This page helps you understand how to use the app through video turotials. Follow the steps below:")
    
    st.subheader("Navigation Tutorial")
    st.video("Videos/navigation.mp4", start_time=0)

    st.subheader("Single Task Tutorial")
    st.video("Videos/start_task1.mp4", start_time=0)
    st.markdown("""
    1. Navigate to the **Start Task** page.
    2. Enter the 8-character code when prompted.
    3. Click the **Start** button to begin the task.
    4. Wait for the task to complete and view the results.""")

    st.subheader("Multi Task Tutorial")
    st.video("Videos/start_task2.mp4", start_time=0)
    st.markdown("""
    1. Navigate to the **Start Task** page.
    2. Enter the 8-character code when prompted.
    3. Click the **Start** button to begin the task.
    4. Wait for the task to complete and view the results.""")

    st.subheader("Glossary Tutorial")
    st.video("Videos/glossary.mp4", start_time=0)

    st.subheader("About Tutorial")
    st.video("Videos/contact us.mp4", start_time=0)
    
# --- About Page ---
elif selected_page == "About":
    st.title("About")
    st.write("**Learn more about the application and its developers.**")

    import urllib.parse
    with st.popover('Contact Us'):
        email_to = "akharbrtk2@gmaill.com"
        subject = "MultiClassClassificationInput App Inquiry"
        body = "Hello, I would like to contact you regarding..."
        subject_encoded = urllib.parse.quote(subject)
        body_encoded = urllib.parse.quote(body)
        # Create the mailto link
        mailto_link = f"mailto:{email_to}?subject={subject_encoded}&body={body_encoded}"
        st.markdown(f"[Tap the link to open E-mail](mailto:{email_to}?subject={subject_encoded}&body={body_encoded})")