INPUT_FILE = "データ/演習/globin.txt"
OUTPUT_FILE = "データ/演習/globin_orf.txt"

PROKARYOTE = 1 << 0
EUKARYOTE = 1 << 1
CATEGORY = PROKARYOTE # 未知配列 => PROKARYOTE | EUKARYOTE

MINIMUM_RESIDUES = 100
CODON_TABLE = [
    "F","F","L","L",
    "S","S","S","S",
    "Y","Y",".",".",
    "C","C",".","W",

    "L","L","L","L",
    "P","P","P","P",
    "H","H","Q","Q",
    "R","R","R","R",

    "I","I","I","M",
    "T","T","T","T",
    "N","N","K","K",
    "S","S","R","R",

    "V","V","V","V",
    "A","A","A","A",
    "D","D","E","E",
    "G","G","G","G"
]

def codon2aa(codon)
    codon_id = 0
    codon.split("").each_with_index do |base,i|
        # T=0,C=1,A=2,G=3とした四進数でコドンを識別
        case base
        when "T"
            codon_id += 0*(4**(2-i))
        when "C"
            codon_id += 1*(4**(2-i))
        when "A"
            codon_id += 2*(4**(2-i))
        when "G"
            codon_id += 3*(4**(2-i))
        end
    end
    return CODON_TABLE[codon_id]
end
def cmpl_strand(strand) # 相補鎖
    strand.split("").map{|base|
        case base
        when "T"
            "A"
        when "C"
            "G"
        when "A"
            "T"
        when "G"
            "C"
        end
    }.join.reverse
end

### 翻訳家 ########
class Translator
    def initialize(seq,category)
        if category & PROKARYOTE == PROKARYOTE
            @results_pro = "-"*25+"{原核生物}"+"-"*25+"\n"
            @orf_counter_pro = 0
        end
        if category & EUKARYOTE == EUKARYOTE
            @results_eu = "-"*25+"{真核生物}"+"-"*25+"\n"
            @orf_counter_eu = 0
        end
        @seq = seq.gsub("\n","") # 配列から改行文字を削除
        @category = category
    end

    def translate_orf(orf,orf_counter,is_cmpl=false)
        results = "【#{orf_counter}】"
        results += "<相補鎖>" if is_cmpl

        ### アミノ酸配列 ########
        results += "\n(アミノ酸配列)\n"
        for i in 0...orf.length/3
            results += codon2aa(orf[3*i,3*(i+1)])+"  " # コドンをアミノ酸に変換
            results += "\n" if i % 20 == 19 # 20残基毎に改行
        end

        orf = cmpl_strand(orf) if is_cmpl # 相補鎖を元の配列に変換

        ### 塩基配列 ########
        results += "\n(塩基配列)\n"
        for i in 0...orf.length
            results += orf[i]
            results += "\n" if i % 60 == 59 # 60塩基毎に改行
        end

        ### 位置 ########
        results += "\n(位置)\n"
        orf_index = @seq.index(orf)
        results += "[#{orf_index},#{orf_index+orf.length-1}]\n\n"
        
        return results
    end

    def translate_all
        [@seq,cmpl_strand(@seq)].each_with_index do |seq_rest,i| # 元の配列と相補鎖
            is_cmpl = (i % 2 == 1)
            while /ATG/ =~ seq_rest # 開始コドン候補を検出
                tmp = $' # 検出された開始コドン以降の配列を保持

                ### 原核生物 ########
                if @category & PROKARYOTE == PROKARYOTE
                    # ORFの候補が存在するか(最短一致)
                    if pseudo_orf = seq_rest.slice(/ATG([TCAG]{3})*?(TAA|TAG|TGA)/)
                        if pseudo_orf.length > 3*MINIMUM_RESIDUES # 配列サイズは十分か
                            @results_pro += translate_orf(pseudo_orf,@orf_counter_pro+=1,is_cmpl)
                        end
                    end
                end

                ### 真核生物 ########
                if @category & EUKARYOTE == EUKARYOTE
                    # ORFの候補が存在するか(最長一致)
                    if pseudo_orf = seq_rest.slice(/ATG([TCAG]{3})*(TAA|TAG|TGA)/)
                        if pseudo_orf.length > 3*MINIMUM_RESIDUES # 配列サイズは十分か
                            @results_eu += translate_orf(pseudo_orf,@orf_counter_eu+=1,is_cmpl)
                        end
                    end
                end

                seq_rest = tmp # 退避変数を代入して反復
            end
        end
        return (@results_pro || "") + (@results_eu || "")
    end
end

### ファイル読み込み ########
file = File.open(INPUT_FILE)
file_content = file.read() # 入力ファイル
new_file_content = file_content # 出力ファイル
file.close

### 翻訳 ########
file_content.scan(/^[TCAG][TCAG\n]+/).each do |seq| # 入力ファイルの塩基配列部分を抽出
    translator = Translator.new(seq,CATEGORY) # 翻訳家に依頼
    new_file_content.sub!(seq,translator.translate_all) # 出力ファイルの塩基配列部分をORFに置換
end

### ファイル書き込み ########
File.open(OUTPUT_FILE,"w") do |text|
    text.puts(new_file_content)
end