#!/bin/bash

# 检查输入参数数量
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <mode> <args...>"
    echo "mode: r (reference-based), d (de novo), or a (all)"
    exit 1
fi

MODE=$1
shift

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# 定义合并GFF文件的函数
merge_gff_files() {
    local gff1="$1"
    local gff2="$2"
    local output_gff="$3"

    # 创建一个临时文件来存储格式化后的第一个 GFF 文件
    temp_gff1="temp_gff1.gff"

    # 1. 格式化第一个 GFF 文件
    while IFS=$'\t' read -r seq source type start end score strand phase attributes; do
    if [[ "$seq" == "#"* ]]; then
        echo "$seq" >> "$temp_gff1"
        continue
    fi

    # 使用 awk 提取属性
    ID=$(echo "$attributes" | awk -F';' '{for(i=1;i<=NF;i++){if($i ~ /^ID=/){print substr($i,4)}}}')
    TSD1_Start=$(echo "$attributes" | awk -F';' '{for(i=1;i<=NF;i++){if($i ~ /^TSD1_Start=/){print substr($i,11)}}}')
    TSD1_End=$(echo "$attributes" | awk -F';' '{for(i=1;i<=NF;i++){if($i ~ /^TSD1_End=/){print substr($i,9)}}}')
    TSD2_Start=$(echo "$attributes" | awk -F';' '{for(i=1;i<=NF;i++){if($i ~ /^TSD2_Start=/){print substr($i,11)}}}')
    TSD2_End=$(echo "$attributes" | awk -F';' '{for(i=1;i<=NF;i++){if($i ~ /^TSD2_End=/){print substr($i,9)}}}')
    TSD1_Sequence=$(echo "$attributes" | awk -F';' '{for(i=1;i<=NF;i++){if($i ~ /^TSD_Sequence=/){print substr($i,13)}}}')
    TSD2_Sequence="$TSD1_Sequence" #  gff1文件中没有TSD2_Sequence，用TSD1_Sequence代替
    TIR1_Start=$(echo "$attributes" | awk -F';' '{for(i=1;i<=NF;i++){if($i ~ /^TIR1_Start=/){print substr($i,11)}}}')
    TIR1_End=$(echo "$attributes" | awk -F';' '{for(i=1;i<=NF;i++){if($i ~ /^TIR1_End=/){print substr($i,9)}}}')
    TIR2_Start=$(echo "$attributes" | awk -F';' '{for(i=1;i<=NF;i++){if($i ~ /^TIR2_Start=/){print substr($i,11)}}}')
    TIR2_End=$(echo "$attributes" | awk -F';' '{for(i=1;i<=NF;i++){if($i ~ /^TIR2_End=/){print substr($i,9)}}}')
    TIR1_Sequence=$(echo "$attributes" | awk -F';' '{for(i=1;i<=NF;i++){if($i ~ /^TIR1_Sequence=/){print substr($i,13)}}}')
    TIR2_Sequence=$(echo "$attributes" | awk -F';' '{for(i=1;i<=NF;i++){if($i ~ /^TIR2_Sequence=/){print substr($i,13)}}}')
    CDS_Start=$(echo "$attributes" | awk -F';' '{for(i=1;i<=NF;i++){if($i ~ /^CDS_Start=/){print substr($i,10)}}}')
    CDS_End=$(echo "$attributes" | awk -F';' '{for(i=1;i<=NF;i++){if($i ~ /^CDS_End=/){print substr($i,8)}}}')


    # 构建新的attributes
    new_attributes="ID=$ID;BLAST_Match=$seq:$start-$end;Element_Structure=$seq:$CDS_Start-$CDS_End;TSD1=$((TSD1_Start - start))-$((TSD1_End - start));TSD1_Sequence=$TSD1_Sequence;TIR1=$((TIR1_Start - start))-$((TIR1_End - start));TIR1_Sequence=$TIR1_Sequence;TSD2=$((TSD2_Start - start))-$((TSD2_End - start));TSD2_Sequence=$TSD2_Sequence;TIR2=$((TIR2_Start - start))-$((TIR2_End - start));TIR2_Sequence=$TIR2_Sequence"

    # 添加 Protein_Region 信息
    if [[ ! -z "$CDS_Start" && ! -z "$CDS_End" ]]; then
        protein_start=$((CDS_Start - start + 1))
        protein_end=$((CDS_End - start))
        new_attributes+=";Protein_Region_1=$protein_start-$protein_end"
    fi

    echo -e "$seq\t$source\t$type\t$start\t$end\t$score\t$strand\t$phase\t$new_attributes" >> "$temp_gff1"

    done < "$gff1"


    # 2. 合并 GFF 文件 (优化版本)

    # 将两个 GFF 文件（包括格式化后的 gff1）合并并排序
    sort -k1,1 -k4,4n "$temp_gff1" "$gff2" > sorted_combined.gff

    # 遍历排序后的合并文件
    prev_seq=""
    prev_start=0
    prev_end=0
    prev_attributes=""
    merged=false

    while IFS=$'\t' read -r seq source type start end score strand phase attributes; do
    if [[ "$seq" == "#"* ]]; then
        echo "$seq" >> "$output_gff"
        continue
    fi

    start_num=$(echo "$start" | awk '{print int($1)}')
    end_num=$(echo "$end" | awk '{print int($1)}')

    if [[ "$seq" == "$prev_seq" ]]; then  # 同一条染色体
        # 计算交集
        intersection_start=$((prev_start > start_num ? prev_start : start_num))
        intersection_end=$((prev_end < end_num ? prev_end : end_num))
        intersection_length=$((intersection_end - intersection_start))

        # 计算并集长度
        union_length=$(( (prev_end > end_num ? prev_end : end_num) - (prev_start < start_num ? prev_start : start_num) ))

        # 计算重叠百分比
        if (( intersection_length > 0 )); then
        overlap_percentage=$(awk -v intersection="$intersection_length" -v union="$union_length" 'BEGIN {printf "%.2f", intersection * 100.0 / union}')
        else
        overlap_percentage=0
        fi

        if (( $(echo "$overlap_percentage > 70" | bc -l) )); then
        # 合并
        merged_start=$((prev_start < start_num ? prev_start : start_num))
        merged_end=$((prev_end > end_num ? prev_end : end_num))
        merged_attributes="$prev_attributes;---merged---;$attributes"
        merged=true
        else
        # 不合并，输出前一个注释
        echo -e "$prev_seq\t$prev_source\t$prev_type\t$prev_start\t$prev_end\t$prev_score\t$prev_strand\t$prev_phase\t$prev_attributes" >> "$output_gff"
        merged=false
        fi
    else  # 不同的染色体，输出前一个注释（如果存在）
        if [[ ! -z "$prev_seq" && ! $merged ]]; then
        echo -e "$prev_seq\t$prev_source\t$prev_type\t$prev_start\t$prev_end\t$prev_score\t$prev_strand\t$prev_phase\t$prev_attributes" >> "$output_gff"
        fi
        merged=false
    fi


    # 更新前一个注释信息
    prev_seq="$seq"
    prev_source="$source"
    prev_type="$type"
    prev_start="$start"
    prev_end="$end"
    prev_score="$score"
    prev_strand="$strand"
    prev_phase="$phase"
    prev_attributes="$attributes"

    done < sorted_combined.gff


    # 输出最后一个注释（如果需要）
    if [[ ! -z "$prev_seq" && ! $merged ]]; then
    echo -e "$prev_seq\t$prev_source\t$prev_type\t$prev_start\t$prev_end\t$prev_score\t$prev_strand\t$prev_phase\t$prev_attributes" >> "$output_gff"
    elif [[ ! -z "$prev_seq" && $merged ]]; then
        echo -e "$prev_seq\t$prev_source\t$prev_type\t$merged_start\t$merged_end\t$prev_score\t$prev_strand\t$prev_phase\t$merged_attributes" >> "$output_gff"
    fi

    rm "$temp_gff1" sorted_combined.gff
    echo "合并完成，输出文件为：$output_gff"
}
if [ "$MODE" = "r" ]; then
    PYTHON_SCRIPT="$SCRIPT_DIR/reference_based_module/reference_based_main.py"
    echo "Running reference-based module..."
    python3 "$PYTHON_SCRIPT" "$@"
    echo "Reference-based module completed."
elif [ "$MODE" = "d" ]; then
    PYTHON_SCRIPT="$SCRIPT_DIR/de_novo_module/transposon_analyzer.py"
    echo "Running de novo module..."
    python3 "$PYTHON_SCRIPT" "$@"
    echo "De novo module completed."
elif [ "$MODE" = "a" ]; then
    echo "Running both reference-based and de novo modules..."

    # 提示用户输入reference-based module需要的参数
    while true; do
        echo "Please enter parameters for the reference-based module, e.g., -g file -o output:"
        read -r ref_based_params
        echo "You entered: $ref_based_params"
        read -p "Is this correct? (yes/no): " confirmation
        if [[ $confirmation =~ ^[Yy] ]]; then
            break
        else
            echo "Please re-enter the parameters."
        fi
    done

    # 提示用户输入de novo module需要的参数
    while true; do
        echo "Please enter parameters for the de novo module, e.g., -i file -o output:"
        read -r de_novo_params
        echo "You entered: $de_novo_params"
        read -p "Is this correct? (yes/no): " confirmation
        if [[ $confirmation =~ ^[Yy] ]]; then
            break
        else
            echo "Please re-enter the parameters."
        fi
    done

    # 提示用户输入最终保存结果的路径
    while true; do
        echo "Please enter the path to save the final merged results, e.g., /path/to/output.gff:"
        read -r output_path
        echo "You will save the results to: $output_path"
        read -p "Is this correct? (yes/no): " confirmation
        if [[ $confirmation =~ ^[Yy] ]]; then
            break
        else
            echo "Please re-enter the output path."
        fi
    done

    # 运行reference-based模块
    echo "Executing reference-based module..."
    if python3 "$SCRIPT_DIR/reference_based_module/reference_based_main.py" $ref_based_params; then
        echo "Reference-based module completed successfully."
    else
        echo "Error running reference-based module. Exiting."
        exit 1
    fi


    echo "Executing de novo module..."
    if python3 "$SCRIPT_DIR/de_novo_module/transposon_analyzer.py" $de_novo_params; then
        echo "De novo module completed successfully."
    else
        echo "Error running de novo module. Exiting."
        exit 1
    fi


    ref_gff_path=$(echo "$ref_based_params" | sed -n 's/.*-gff \(\S\+\).*/\1/p')
    de_novo_gff_path=$(echo "$de_novo_params" | sed -n 's/.*-o \(\S\+\).*/\1/p')


    echo "Merging results from both modules..."
    merge_gff_files "$ref_gff_path" "$de_novo_gff_path" "$output_path"
    echo "All modules completed and results merged."
else
    echo "Invalid mode. Use 'r' for reference-based, 'd' for de novo, or 'a' for all."
    exit 1
fi