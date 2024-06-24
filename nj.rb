require "set"
require "csv"
require "bigdecimal"
flag = BigDecimal::ROUND_DOWN
f = BigDecimal::mode(BigDecimal::ROUND_MODE,flag) # 丸め処理を切り捨てに指定

INPUT_FILE = "データ/設問/aldolase2_cdis.csv"
OUTPUT_FILE = "データ/設問/aldolase2_NJ.txt"

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

### 近隣結合法 ########
while entries.length > 2
    n = entries.length # OTUの数
    s_min = Float::INFINITY # 枝の総和の最小値
    min_pair = [] ; min_pair_l = [] # 最小ペアと枝の長さ

    ### 最小ペアの検出 ########
    cdis.each_key do |name_pair|
        i,j = name_pair.split(" & ") # 仮結合ペア
        r_i = entries.dup.delete(i).map{|entry| cdis[i+" & "+entry] || cdis[entry+" & "+i]}.sum
        r_j = entries.dup.delete(j).map{|entry| cdis[j+" & "+entry] || cdis[entry+" & "+j]}.sum
        s_ij = (n-2)*cdis[name_pair] - r_i - r_j # 枝の総和の簡易値
        if s_ij < s_min
            s_min = s_ij # 最小値の更新
            min_pair = [i,j]
            # 負の長さ => 0
            l_i = [round(((n-2)*cdis[name_pair] + r_i - r_j)/(2*(n-2))),0].max
            l_j = [round(((n-2)*cdis[name_pair] - r_i + r_j)/(2*(n-2))),0].max
            min_pair_l = [l_i,l_j]
        end
    end
    cdis_min = [min_pair.join(" & "),min_pair_l]
    results.push(cdis_min) # 最小ペア名と枝の長さを記録

    ### 距離行列と集合の更新 ########
    entries.delete(min_pair[0]) ; entries.delete(min_pair[1]) # 集合から各配列を削除
    entries.each do |entry|
        new_key = "(#{min_pair[0]} + #{min_pair[1]}) & #{entry}" # 最小ペアを結合ノードに変換
        # 結合前ペア
        key1_1 = min_pair[0]+" & "+entry ; key1_2 = entry+" & "+min_pair[0]
        key2_1 = min_pair[1]+" & "+entry ; key2_2 = entry+" & "+min_pair[1]
        value1 = cdis[key1_1] || cdis[key1_2] ; value2 = cdis[key2_1] || cdis[key2_2]
        # 結合ノードで連想配列を更新
        cdis[new_key] = [round((value1 + value2 - cdis[results.last[0]])/2.0),0].max # 負の距離 => 0
        # 結合前ペアを削除
        cdis.delete(key1_1) ; cdis.delete(key1_2) ; cdis.delete(key2_1) ; cdis.delete(key2_2)
    end
    cdis.delete(cdis_min[0]) # 最小ペアを連想配列から削除
    entries.add("(#{min_pair[0]} + #{min_pair[1]})") # 集合に結合ノードを追加

end
results.push(cdis) # 残ノード間距離を記録

### ファイル書き込み ########
File.open(OUTPUT_FILE,"w") do |text|
    results.each_with_index do |line,i|
        text.puts("STEP#{i} #{line}")
    end
end