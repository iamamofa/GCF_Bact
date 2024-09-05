import customtkinter as ctk
from tkinter import filedialog, messagebox
from pathlib import Path
import subprocess

ctk.set_appearance_mode("Dark")     # Modes: "System", "Dark", "Light"
ctk.set_default_color_theme("green")    # Themes: "blue", "dark-blue", "green"

class App(ctk.CTk):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.title("Bacterial Genomics Pipeline")
        self.geometry("1020x1080")
        self.configure(bg="#1f1f1f")

        # Create a sidebar frame
        self.sidebar_frame = ctk.CTkFrame(self, width=200, corner_radius=10)
        self.sidebar_frame.grid(row=0, column=0, rowspan=4, sticky="nswe", padx=20, pady=20)

        # Sidebar label
        self.sidebar_label = ctk.CTkLabel(self.sidebar_frame, text="Pipeline Options", font=("Arial", 20))
        self.sidebar_label.grid(row=0, column=0, padx=20, pady=(20, 10))

        # Create directory and file selection widgets
        self.directory_button = ctk.CTkButton(self.sidebar_frame, text="Select Result Directory", command=self.select_directory)
        self.directory_button.grid(row=1, column=0, padx=20, pady=10)

        self.file_button = ctk.CTkButton(self.sidebar_frame, text="Upload FASTA/FASTQ File", command=self.upload_file)
        self.file_button.grid(row=2, column=0, padx=20, pady=10)

        self.run_button = ctk.CTkButton(self.sidebar_frame, text="Run Pipeline", command=self.run_pipeline)
        self.run_button.grid(row=3, column=0, padx=20, pady=(10, 20))

        # Create a main content frame
        self.main_frame = ctk.CTkFrame(self, corner_radius=10)
        self.main_frame.grid(row=0, column=1, rowspan=4, sticky="nswe", padx=20, pady=20)

        # Progress section
        self.progress_label = ctk.CTkLabel(self.main_frame, text="Pipeline Progress", font=("Arial", 18))
        self.progress_label.pack(pady=20)
        
        self.progress_text = ctk.CTkTextbox(self.main_frame, width=760, height=300)
        self.progress_text.pack(pady=20)
        
        # Result directory and uploaded file placeholders
        self.result_dir = None
        self.uploaded_file = None

    def select_directory(self):
        self.result_dir = Path(filedialog.askdirectory())
        if self.result_dir:
            self.progress_text.insert(ctk.END, f"Selected directory: {self.result_dir}\n")
            messagebox.showinfo("Selected Directory", str(self.result_dir))

    def upload_file(self):
        filetypes = [('FASTA/FASTQ Files', '*.fasta *.fastq')]
        self.uploaded_file = Path(filedialog.askopenfilename(filetypes=filetypes))
        if self.uploaded_file:
            self.progress_text.insert(ctk.END, f"Uploaded file: {self.uploaded_file}\n")
            messagebox.showinfo("Uploaded File", str(self.uploaded_file))

    def run_pipeline(self):
        if not self.result_dir or not self.uploaded_file:
            messagebox.showwarning("Input Error", "Please select a directory and upload a file.")
            return
        
        # Determine file type
        file_extension = self.uploaded_file.suffix.lower()
        
        if file_extension == ".fasta":
            self.run_index_reference()
        elif file_extension == ".fastq":
            self.run_quality_adapter_trimming()
        else:
            messagebox.showerror("File Error", "Invalid file type. Please upload a FASTA or FASTQ file.")

    def run_index_reference(self):
        self.progress_text.insert(ctk.END, "Running reference indexing...\n")
        command = ["bash", "./src/index_reference.sh", str(self.uploaded_file), str(self.result_dir)]
        self._run_command(command)
        self.progress_text.insert(ctk.END, "Reference indexing complete.\n")
        self.run_map_reads()

    def run_quality_adapter_trimming(self):
        self.progress_text.insert(ctk.END, "Running quality and adapter trimming...\n")
        command = ["bash", "./src/quality_adapter_trimming.sh", str(self.uploaded_file), str(self.result_dir)]
        self._run_command(command)
        self.progress_text.insert(ctk.END, "Trimming complete.\n")
        self.run_map_reads()

    def run_map_reads(self):
        self.progress_text.insert(ctk.END, "Mapping reads...\n")
        indexed_reference = self.result_dir / "indexed_reference"
        command = ["bash", "./src/map_reads.sh", str(self.uploaded_file), str(indexed_reference), str(self.result_dir)]
        self._run_command(command)
        self.progress_text.insert(ctk.END, "Reads mapping complete.\n")
        self.run_sort_index_reads()

    def run_sort_index_reads(self):
        self.progress_text.insert(ctk.END, "Sorting and indexing reads...\n")
        input_sam = self.result_dir / "aligned_reads.sam"
        command = ["bash", "./src/sort_index_reads.sh", str(input_sam), str(self.result_dir)]
        self._run_command(command)
        self.progress_text.insert(ctk.END, "Sorting and indexing complete.\n")
        self.run_variant_calling()

    def run_variant_calling(self):
        self.progress_text.insert(ctk.END, "Running variant calling...\n")
        sorted_bam = self.result_dir / "sorted_reads.bam"
        reference = self.result_dir / "indexed_reference"
        command = ["bash", "./src/variant_calling.sh", str(sorted_bam), str(reference), str(self.result_dir)]
        self._run_command(command)
        self.progress_text.insert(ctk.END, "Variant calling complete.\n")
        self.run_filter_variants()

    def run_filter_variants(self):
        self.progress_text.insert(ctk.END, "Filtering variants...\n")
        input_vcf = self.result_dir / "variants_bcftools.vcf.gz"
        command = ["bash", "./src/filter_variants.sh", str(input_vcf), str(self.result_dir)]
        self._run_command(command)
        self.progress_text.insert(ctk.END, "Filtering complete.\n")
        self.run_create_consensus_fasta()

    def run_create_consensus_fasta(self):
        self.progress_text.insert(ctk.END, "Creating consensus FASTA...\n")
        input_vcf = self.result_dir / "filtered_variants.vcf"
        reference = self.result_dir / "indexed_reference"
        command = ["bash", "./src/create_consensus_fasta.sh", str(input_vcf), str(reference), str(self.result_dir)]
        self._run_command(command)
        self.progress_text.insert(ctk.END, "Consensus FASTA creation complete.\n")
        self.run_remove_non_informative_sites()

    def run_remove_non_informative_sites(self):
        self.progress_text.insert(ctk.END, "Removing non-informative sites...\n")
        input_fasta = self.result_dir / "consensus.fasta"
        command = ["bash", "./src/remove_non_informative_sites.sh", str(input_fasta), str(self.result_dir)]
        self._run_command(command)
        self.progress_text.insert(ctk.END, "Non-informative sites removal complete.\n")
        self.run_build_tree()

    def run_build_tree(self):
        self.progress_text.insert(ctk.END, "Building phylogenetic tree...\n")
        input_fasta = self.result_dir / "informative_sites.fasta"
        command = ["bash", "./src/build_tree.sh", str(input_fasta), str(self.result_dir)]
        self._run_command(command)
        self.progress_text.insert(ctk.END, "Tree building complete. Check the results directory for outputs.\n")
        messagebox.showinfo("Pipeline Complete", "All steps completed successfully!")

    def _run_command(self, command):
        try:
            result = subprocess.run(command, check=True, text=True, capture_output=True)
            self.progress_text.insert(ctk.END, f"Command output:\n{result.stdout}\n")
        except subprocess.CalledProcessError as e:
            self.progress_text.insert(ctk.END, f"Error running command:\n{e.stderr}\n")
            messagebox.showerror("Command Error", f"An error occurred: {e.stderr}")

if __name__ == "__main__":
    app = App()
    app.mainloop()
