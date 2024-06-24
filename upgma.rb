require "set"
require "csv"
require "bigdecimal"
flag = BigDecimal::ROUND_DOWN
f = BigDecimal::mode(BigDecimal::ROUND_MODE,flag) # 丸め処理を切り捨てに指定

INPUT_FILE = "データ/設問/aldolase2_cdis.csv"
OUTPUT_FILE = "データ/設問/aldolase2_UPGMA.txt"

def round(value)
    BigDecimal("0").add(value,4).to_f # 有効数字3桁+1桁
end

### ファイル読み込み ########
cdis = {} # 配列ペアの進化距離の連想配列(距離行列)
entries = Set.new # 配列名の集合
CSV.foreach(INPUT_FILE) do |row|
    cdis[row[0].strip] = row[1].to_f
    name_pair = row[0].strip.split(" & ")
    entries.add(name_pair[0]) ; entries.add(name_pair[1])
end
results = []

### 平均距離法 ########
while entries.length > 1

    ### 最小ペアの検出 ########
    cdis_min = cdis.min{|x,y| x[1] <=> y[1]}
    results.push([cdis_min[0],round(cdis_min[1]/2.0)]) # 最小ペア名と枝の長さ(=進化距離/2)を記録

    ### 距離行列と集合の更新 ########
    min_pair = cdis_min[0].split(" & ") # 最小ペアを分割
    entries.delete(min_pair[0]) ; entries.delete(min_pair[1]) # 集合から各配列を削除
    entries.each do |entry|
        new_key = "(#{min_pair[0]} + #{min_pair[1]}) & #{entry}" # 最小ペアを結合ノードに変換
        # 結合前ペア
        key1_1 = min_pair[0]+" & "+entry ; key1_2 = entry+" & "+min_pair[0]
        key2_1 = min_pair[1]+" & "+entry ; key2_2 = entry+" & "+min_pair[1]
        value1 = cdis[key1_1] || cdis[key1_2] ; value2 = cdis[key2_1] || cdis[key2_2]
        # 結合ノードで連想配列を更新
        cdis[new_key] = round((value1 + value2)/2.0)
        # 結合前ペアを削除
        cdis.delete(key1_1) ; cdis.delete(key1_2) ; cdis.delete(key2_1) ; cdis.delete(key2_2)
    end
    cdis.delete(cdis_min[0]) # 最小ペアを連想配列から削除
    entries.add("(#{min_pair[0]} + #{min_pair[1]})") # 集合に結合ノードを追加
    
end

### ファイル書き込み ########
File.open(OUTPUT_FILE,"w") do |text|
    results.each_with_index do |line,i|
        text.puts("STEP#{i} #{line}")
    end
end