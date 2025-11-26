current_dir=$(pwd)
awk '$5<0.2{
    system("cd "$1" && \
    if [ -e sp ]; then \
        echo \"Skipping "$1": sp file already exists\"; \
    else \
        echo \"Processing "$1"\"; \
        align_lines=$(grep -c \"^\" ALIGN 2>/dev/null || echo 0); \
        n_value=$(( align_lines < 100 ? align_lines : 100 )); \
        if [ $n_value -gt 0 ]; then \
            swissprot.pl -n $n_value; \
        else \
            echo \"Warning: No ALIGN file or empty ALIGN in "$1"\"; \
        fi; \
    fi && \
    cd '"$current_dir"'")
}' triage
