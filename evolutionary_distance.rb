require "bigdecimal"
flag = BigDecimal::ROUND_DOWN
f = BigDecimal::mode(BigDecimal::ROUND_MODE,flag) # 丸め処理を切り捨てに指定

INPUT_FILE = "データ/設問/aldolase2.txt"
OUTPUT_FILE = "データ/設問/aldolase2_cdis.csv"

### ファイル読み込み ########
file = File.open(INPUT_FILE)
file_content = file.read() # 入力ファイル
seq_names = file_content.scan(/^=.*/).map{|header| header.delete("=").strip} # ヘッダー部分の抽出と加工
seq_name_pairs = seq_names.combination(2).to_a # 配列名の組を生成
seqs = file_content.scan(/^[[A-Z]-][[A-Z]-\n]+/).map{|seq| seq.gsub("\n","")} # 配列部分の抽出と加工
seq_pairs = seqs.combination(2).to_a # 配列の組を生成
ict = 0 # 有効座位数
sbst = Array.new(seq_pairs.length,0) # 各組の置換数の配列
file.close

### 進化距離計算 ########
for i in 0...seqs[0].length # 第i座位について
    next if i == 0 && seqs[0][i] == "M" # 第0座位が開始メチオニンならスキップ
    next if seqs.any?{|seq| seq[i] == "-"} # ギャップが1つでもあればスキップ
    ict += 1
    seq_pairs.each_with_index do |seq_pair,j|
        sbst[j] += 1 if seq_pair[0][i] != seq_pair[1][i] # 配列間でアミノ酸が異なれば置換数を増分
    end
end
# 相違度のPoisson補正値(進化距離)の配列
cdis = sbst.map{|s| BigDecimal("0").add(-Math.log(1.0 - s/ict.to_f),4).to_f} # 有効数字3桁+1桁

### ファイル書き込み ########
File.open(OUTPUT_FILE,"w") do |text|
    for k in 0...seq_pairs.length
        text.puts(seq_name_pairs[k][0]+" & "+seq_name_pairs[k][1]+" , "+cdis[k].to_s)
    end
end