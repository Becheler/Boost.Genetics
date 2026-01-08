<invoke name="artifacts">
<parameter name="command">update</parameter>
<parameter name="id">vcf_impl_spec</parameter>
<parameter name="old_str">### Example 2.4: Insertion

```vcf
##</parameter>
<parameter name="new_str">### Example 2.4: Insertion

```vcf
##fileformat=VCFv4.3
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
chr1	400	.	C	CTAG	50	PASS	.	GT	0/1
```

**Expected Parsing:**
```
Variant Type: Insertion
Inserted Bases: "TAG" (after position 400)
Padding Base: "C" (at position 400)
Sample1 Genotype: Heterozygous insertion
```

---

### Example 2.5: Missing Genotype Values

```vcf
##fileformat=VCFv4.3
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	500	.	A	G	30	PASS	.	GT:GQ:DP	0/1:.:10	./.:20:.
```

**Expected Parsing:**
```
Samples:
  S1:
    GT: "0/1"
    GQ: MISSING
    DP: 10
  S2:
    GT: "./." (both alleles missing)
    GQ: 20
    DP: MISSING
```

---

### Example 2.6: Phased Genotypes

```vcf
##fileformat=VCFv4.3
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
chr1	100	.	A	G	30	PASS	.	GT:PS	0|1:100
chr1	200	.	C	T	30	PASS	.	GT:PS	1|0:100
chr1	300	.	G	A	30	PASS	.	GT	0/1
```

**Expected Parsing:**
```
Record 1 & 2: Same phase set (PS=100), phased together
Record 3: Unphased (no PS value)
```

---

### Example 2.7: Variable-Length Genotype Fields

```vcf
##fileformat=VCFv4.3
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled likelihoods">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	600	.	A	G,T	50	PASS	.	GT:PL	0/1:10,0,20,.,.,.	1/2:.,.,.,5,0,10
```

**Expected Parsing:**
```
N_ALT = 2, so G = 6 genotypes (0/0, 0/1, 1/1, 0/2, 1/2, 2/2)

S1:
  GT: "0/1"
  PL: [10, 0, 20, MISSING, MISSING, MISSING]
  # Only first 3 values provided

S2:
  GT: "1/2"
  PL: [MISSING, MISSING, MISSING, 5, 0, 10]
  # Only last 3 values provided
```

---

## Phase 3: Structural Variants ✅

### 3.1 Reserved INFO Keys for SVs

**Spec Reference:** Section 3

| Key       | Number | Type    | Description |
|-----------|--------|---------|-------------|
| IMPRECISE | 0      | Flag    | Imprecise structural variation |
| NOVEL     | 0      | Flag    | Novel structural variation |
| END       | 1      | Integer | End position of variant |
| SVTYPE    | 1      | String  | Type of SV (DEL/INS/DUP/INV/CNV/BND) |
| SVLEN     | .      | Integer | Difference in length (REF - ALT) |
| CIPOS     | 2      | Integer | Confidence interval around POS |
| CIEND     | 2      | Integer | Confidence interval around END |
| HOMLEN    | .      | Integer | Length of micro-homology |
| HOMSEQ    | .      | String  | Sequence of micro-homology |
| BKPTID    | .      | String  | ID in assembly file |
| MEINFO    | 4      | String  | Mobile element info |
| METRANS   | 4      | String  | Mobile element transduction |
| DGVID     | 1      | String  | Database of Genomic Variation ID |
| DBVARID   | 1      | String  | DBVAR ID |
| MATEID    | .      | String  | ID of mate breakend |
| PARID     | 1      | String  | ID of partner breakend |
| EVENT     | 1      | String  | Event ID |
| CILEN     | 2      | Integer | CI around inserted material |
| DPADJ     | .      | Integer | Read depth of adjacency |
| CNADJ     | .      | Integer | Copy number of adjacency |

---

### 3.2 Reserved FORMAT Keys for SVs

**Spec Reference:** Section 4

| Key  | Number | Type    | Description |
|------|--------|---------|-------------|
| CN   | 1      | Integer | Copy number genotype |
| CNQ  | 1      | Float   | Copy number quality |
| CNL  | G      | Float   | Copy number likelihoods |
| CNP  | G      | Float   | Copy number posterior probabilities |
| NQ   | 1      | Integer | Phred probability variant is novel |
| HAP  | 1      | Integer | Haplotype identifier |
| AHAP | 1      | Integer | Ancestral haplotype identifier |

---

### 3.3 Simple Deletion Example

```vcf
##fileformat=VCFv4.3
##ALT=<ID=DEL,Description="Deletion">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length difference">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
chr1	1000	.	N	<DEL>	.	PASS	SVTYPE=DEL;END=2000;SVLEN=-1000	GT	1/1
```

**Expected Parsing:**
```
Variant Type: Deletion
Start: 1000
End: 2000
Length: 1000 bp deleted
Genotype: Homozygous deletion
```

---

### 3.4 Imprecise Deletion

```vcf
##fileformat=VCFv4.3
##ALT=<ID=DEL,Description="Deletion">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length difference">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
chr1	5000	.	N	<DEL>	6	PASS	IMPRECISE;SVTYPE=DEL;END=5205;SVLEN=-205;CIPOS=-56,20;CIEND=-10,62	GT	0/1
```

**Expected Parsing:**
```
Variant Type: Imprecise Deletion
Approximate Start: 5000 (range: 4944-5020)
Approximate End: 5205 (range: 5195-5267)
Approximate Length: ~205 bp
Genotype: Heterozygous
```

---

### 3.5 Breakend Notation

**Spec Reference:** Section 5.4

**Breakend Replacement Formats:**

| REF | ALT       | Meaning |
|-----|-----------|---------|
| s   | t[p[      | Piece extending right of p is joined after t |
| s   | t]p]      | Reverse complement piece extending left of p is joined after t |
| s   | ]p]t      | Piece extending left of p is joined before t |
| s   | [p[t      | Reverse complement piece extending right of p is joined before t |

**Where:**
- `s` = bases being replaced
- `t` = replacement bases (possibly with novel inserted sequence)
- `p` = mate breakend position (format: "chr:pos")

**Example:**
```vcf
##fileformat=VCFv4.3
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">
##INFO=<ID=MATEID,Number=.,Type=String,Description="Mate breakend ID">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr2	321682	bnd_V	T	]chr13:123456]T	6	PASS	SVTYPE=BND;MATEID=bnd_U
chr13	123456	bnd_U	C	C[chr2:321682[	6	PASS	SVTYPE=BND;MATEID=bnd_V
```

**Expected Parsing:**
```
Two breakends forming a novel adjacency:
  bnd_V at chr2:321682 joins to bnd_U at chr13:123456
```

---

## Phase 4: gVCF Support ✅

### 4.1 Reference Blocks

**Spec Reference:** Section 5.5

**Purpose:** Represent evidence for non-variant positions

**Format:** Use `<*>` symbolic allele for "any other allele"

**Example:**
```vcf
##fileformat=VCFv4.3
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP in block">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled likelihoods">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
chr1	100	.	G	<*>	.	.	END=110	GT:DP:GQ:MIN_DP:PL	0/0:25:60:23:0,60,900
chr1	111	.	T	C,<*>	213	.	.	GT:DP:GQ:PL	0/1:23:99:51,0,36,93,92,86
chr1	112	.	C	<*>	.	.	END=120	GT:DP:GQ:MIN_DP:PL	0/0:27:63:27:0,63,945
```

**Expected Parsing:**
```
Positions 100-110: Reference blocks (no variant)
  - Minimum depth in block: 23
  - Quality: 60

Position 111: Variant (T→C)
  - Heterozygous

Positions 112-120: Reference blocks (no variant)
  - Minimum depth in block: 27
  - Quality: 63
```

---

## Phase 5: BCF Binary Format (Advanced)

### 5.1 BCF File Structure

**Spec Reference:** Section 6

**Components:**
1. BCF Header
2. BGZF-compressed BCF Records
3. Empty BGZF block (EOF marker)

---

### 5.2 BCF Header Format

**Spec Reference:** Section 6.2

```
Offset | Type     | Field      | Description
-------|----------|------------|-------------
0      | char[3]  | magic      | "BCF"
3      | uint8_t  | major_ver  | 2
4      | uint8_t  | minor_ver  | 2
5      | uint32_t | l_text     | Length of VCF header text
9      | char[]   | text       | VCF header (NUL-terminated)
```

**Test Case 5.2.1:**
```cpp
// Read BCF header
assert(header[0] == 'B');
assert(header[1] == 'C');
assert(header[2] == 'F');
assert(header[3] == 2);  // major version
assert(header[4] == 2);  // minor version
```

---

### 5.3 BCF Type Encoding

**Spec Reference:** Section 6.3.3

**Type Descriptor Byte:**
```
Bits 5-8: Number of elements (or 15 for overflow)
Bits 1-4: Type code
```

**Type Codes:**
```
0x0 = MISSING
0x1 = INT8
0x2 = INT16
0x3 = INT32
0x5 = FLOAT
0x7 = CHAR
```

**Reserved Integer Values:**
```
INT8:  0x80 = MISSING, 0x81 = END_OF_VECTOR, 0x82-0x87 = Reserved
INT16: 0x8000 = MISSING, 0x8001 = END_OF_VECTOR, 0x8002-0x8007 = Reserved
INT32: 0x80000000 = MISSING, 0x80000001 = END_OF_VECTOR, 0x80000002-0x80000007 = Reserved
```

**Reserved Float Values:**
```
0x7F800001 = MISSING
0x7F800002 = END_OF_VECTOR
0x7F800003-0x7F800007 = Reserved
0x7FC00000 = NaN (quiet NaN)
```

**Examples:**

**Atomic INT8 with value 42:**
```
0x11 0x2A
```

**3-element INT8 vector [1, 2, 3]:**
```
0x31 0x01 0x02 0x03
```

**String "ACGT":**
```
0x47 0x41 0x43 0x47 0x54
(0x40 | 0x07 = type, then ASCII bytes)
```

**Vector with overflow size (16 elements):**
```
0xF1           # Type with size=15 (overflow)
0x11 0x10      # Actual size as typed INT8 (16)
...16 bytes... # Actual data
```

---

### 5.4 BCF GT Encoding

**Spec Reference:** Section 6.4.7

**Formula:** `(allele + 1) << 1 | phased`

**Examples:**

| VCF GT | Binary Encoding | Hex     |
|--------|-----------------|---------|
| 0/0    | [0x02, 0x02]    | 02 02   |
| 0/1    | [0x02, 0x04]    | 02 04   |
| 1/1    | [0x04, 0x04]    | 04 04   |
| 0\|1   | [0x02, 0x05]    | 02 05   |
| 1\|0   | [0x05, 0x02]    | 05 02   |
| ./.    | [0x00, 0x00]    | 00 00   |
| 0      | [0x02]          | 02      |
| 1      | [0x04]          | 04      |
| 0/1/2  | [0x02,0x04,0x06]| 02 04 06|

**Mixed Ploidy Example:**
```
Sample1: "0" (haploid)
Sample2: "0/1" (diploid)

Encoded:
Sample1: [0x02, 0x81]  # 0x81 = END_OF_VECTOR
Sample2: [0x02, 0x04]
```

---

## Phase 6: Validation Rules

### 6.1 Header Validation

**Rules:**
1. `##fileformat=VCFv4.3` MUST be first line
2. All `##INFO`, `##FORMAT`, `##FILTER` lines must have unique IDs
3. All IDs must match: `^([A-Za-z_][0-9A-Za-z_.]*)$` (except "1000G")
4. No duplicate sample names in header
5. FORMAT column required if samples present
6. All contigs referenced in data must be defined in header

**Test Cases:**

**Test 6.1.1:** (Missing fileformat)
```vcf
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
```
Expected: ERROR

**Test 6.1.2:** (Duplicate INFO ID)
```vcf
##fileformat=VCFv4.3
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=DP,Number=1,Type=Float,Description="Different Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
```
Expected: ERROR

---

### 6.2 Data Line Validation

**Rules:**
1. Must have exactly 8 or 8+N columns (N = FORMAT + samples)
2. POS must be positive integer (or 0 for telomere)
3. REF cannot be empty
4. ALT must match: `^([ACGTNacgtn]+|\*|\.)$` or be symbolic `<ID>`
5. All FILTER values must be defined in header (except "PASS" and ".")
6. All INFO keys must be defined in header
7. INFO values must match declared types
8. All FORMAT keys must be defined in header
9. Number of sample fields must match FORMAT field count (or fewer if trailing)

**Test Cases:**

**Test 6.2.1:** (Invalid ALT)
```vcf
##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	X	30	PASS	.
```
Expected: ERROR (X not valid base)

**Test 6.2.2:** (Undefined INFO key)
```vcf
##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	G	30	PASS	UNDEFINED=50
```
Expected: ERROR (UNDEFINED not in header)

**Test 6.2.3:** (Type mismatch)
```vcf
##fileformat=VCFv4.3
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	G	30	PASS	DP=abc
```
Expected: ERROR ("abc" not an integer)

---

### 6.3 Sorting Validation

**Rules:**
1. Records for same CHROM must be contiguous
2. Within each CHROM, POS must be non-decreasing
3. Multiple records at same POS are allowed

**Test Case 6.3.1:** (Out of order)
```vcf
##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	200	.	A	G	30	PASS	.
chr1	100	.	C	T	30	PASS	.
```
Expected: ERROR (POS not sorted)

**Test Case 6.3.2:** (Non-contiguous CHROM)
```vcf
##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	G	30	PASS	.
chr2	200	.	C	T	30	PASS	.
chr1	300	.	G	A	30	PASS	.
```
Expected: ERROR (chr1 not contiguous)

---

## Phase 7: Performance Considerations

### 7.1 Parsing Optimizations

**Tips:**
1. Use string views to avoid allocations
2. Parse numbers directly without string conversion when possible
3. Reserve space for vectors based on Number field
4. Cache header dictionaries for O(1) lookup
5. Use memory-mapped I/O for large files

**Example:**
```cpp
// Bad: Multiple allocations
std::vector<std::string> split(const std::string& str, char delim) {
    std::vector<std::string> result;
    std::stringstream ss(str);
    std::string item;
    while (std::getline(ss, item, delim)) {
        result.push_back(item);
    }
    return result;
}

// Good: Use string_view
std::vector<std::string_view> split_view(std::string_view str, char delim) {
    std::vector<std::string_view> result;
    size_t start = 0;
    size_t end = str.find(delim);
    while (end != std::string_view::npos) {
        result.push_back(str.substr(start, end - start));
        start = end + 1;
        end = str.find(delim, start);
    }
    result.push_back(str.substr(start));
    return result;
}
```

---

### 7.2 Memory Management

**Strategies:**
1. **Streaming:** Process one record at a time
2. **Batching:** Load records in chunks
3. **Indexing:** Use tabix/CSI for random access
4. **Lazy Parsing:** Only parse fields when accessed

---

## Phase 8: Testing Strategy

### 8.1 Unit Tests

**Test Categories:**
1. **Type Parsing:** All data types with edge cases
2. **Encoding/Decoding:** Percent encoding, GT encoding
3. **Header Parsing:** All meta-information types
4. **Data Parsing:** Fixed fields, INFO, FORMAT
5. **Validation:** All error conditions
6. **Edge Cases:** Empty fields, missing values, large values

---

### 8.2 Integration Tests

**Test Files:**
1. **Minimal:** Smallest valid VCF
2. **Standard:** Common variants (SNPs, indels)
3. **Complex:** Multi-allelic, phased, structural variants
4. **gVCF:** Reference blocks
5. **Invalid:** Various error conditions
6. **Real Data:** 1000 Genomes, gnomAD samples

---

### 8.3 Comprehensive Test Suite

**Test File 8.3.1: Minimal Valid**
```vcf
##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	1	.	A	T	.	.	.
```

**Test File 8.3.2: All Field Types**
```vcf
##fileformat=VCFv4.3
##INFO=<ID=INT,Number=1,Type=Integer,Description="Integer">
##INFO=<ID=FLOAT,Number=1,Type=Float,Description="Float">
##INFO=<ID=STRING,Number=1,Type=String,Description="String">
##INFO=<ID=CHAR,Number=1,Type=Character,Description="Character">
##INFO=<ID=FLAG,Number=0,Type=Flag,Description="Flag">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	1	.	A	T	.	.	INT=42;FLOAT=3.14;STRING=test;CHAR=X;FLAG
```

**Test File 8.3.3: All Special Values**
```vcf
##fileformat=VCFv4.3
##INFO=<ID=TEST,Number=1,Type=Float,Description="Test">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	1	.	A	T	.	.	TEST=.
1	2	.	A	T	INF	.	TEST=INF
1	3	.	A	T	-INF	.	TEST=-INF
1	4	.	A	T	NAN	.	TEST=NAN
```

**Test File 8.3.4: Complex Genotypes**
```vcf
##fileformat=VCFv4.3
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled likelihoods">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
1	1	.	A	G,T	.	.	.	GT:PL	0/0:0,10,100,10,100,100	0/1:10,0,100,10,100,100	1/2:100,100,100,10,0,10	2/2:100,100,100,100,100,0
```

---

## Phase 9: API Design

### 9.1 Core Classes

```cpp
// Pseudo-code for header-only library structure

class VCFHeader {
public:
    struct InfoField {
        std::string id;
        std::string number;  // "1", "A", "R", "G", "."
        std::string type;    // "Integer", "Float", "String", "Character", "Flag"
        std::string description;
        std::optional<std::string> source;
        std::optional<std::string> version;
    };
    
    struct FormatField {
        std::string id;
        std::string number;
        std::string type;
        std::string description;
    };
    
    struct FilterField {
        std::string id;
        std::string description;
    };
    
    struct ContigField {
        std::string id;
        std::optional<int> length;
        std::optional<std::string> md5;
        std::optional<std::string> url;
    };
    
    int major_version;
    int minor_version;
    std::vector<std::string> samples;
    std::unordered_map<std::string, InfoField> info_fields;
    std::unordered_map<std::string, FormatField> format_fields;
    std::unordered_map<std::string, FilterField> filter_fields;
    std::unordered_map<std::string, ContigField> contigs;
    std::vector<std::string> extra_lines;  // Other meta-information
};

class VCFRecord {
public:
    struct InfoValue {
        std::string key;
        std::variant<int, float, std::string, bool, std::vector<int>, std::vector<float>, std::vector<std::string>> value;
    };
    
    struct GenotypeValue {
        std::string key;
        std::variant<int, float, std::string, std::vector<int>, std::vector<float>, std::vector<std::string>> value;
    };
    
    struct Sample {
        std::string name;
        std::vector<GenotypeValue> values;
        
        std::optional<std::string> get_gt() const;
        std::optional<int> get_gq() const;
        std::optional<int> get_dp() const;
        // ... other helper methods
    };
    
    std::string chrom;
    int pos;
    std::vector<std::string> id;
    std::string ref;
    std::vector<std::string> alt;
    std::optional<float> qual;
    std::vector<std::string> filter;
    std::vector<InfoValue> info;
    std::vector<Sample> samples;
    
    // Helper methods
    bool is_snp() const;
    bool is_indel() const;
    bool is_deletion() const;
    bool is_insertion() const;
    bool is_sv() const;
    int variant_length() const;
};

class VCFParser {
public:
    VCFParser(std::istream& input);
    
    const VCFHeader& header() const;
    bool next(VCFRecord& record);
    
private:
    std::istream& input_;
    VCFHeader header_;
    std::string line_buffer_;
};

class VCFWriter {
public:
    VCFWriter(std::ostream& output, const VCFHeader& header);
    
    void write(const VCFRecord& record);
    
private:
    std::ostream& output_;
    VCFHeader header_;
};
```

---

### 9.2 Usage Examples

**Example 9.2.1: Basic Parsing**
```cpp
std::ifstream file("variants.vcf");
VCFParser parser(file);

const auto& header = parser.header();
std::cout << "Samples: ";
for (const auto& sample : header.samples) {
    std::cout << sample << " ";
}
std::cout << "\n";

VCFRecord record;
while (parser.next(record)) {
    std::cout << record.chrom << ":" << record.pos << " "
              << record.ref << "->" << record.alt[0] << "\n";
}
```

**Example 9.2.2: Filtering Variants**
```cpp
VCFParser parser(std::cin);
VCFWriter writer(std::cout, parser.header());

VCFRecord record;
while (parser.next(record)) {
    // Only output high-quality SNPs
    if (record.is_snp() && record.qual && *record.qual >= 30) {
        writer.write(record);
    }
}
```

**Example 9.2.3: Accessing Genotypes**
```cpp
VCFParser parser(file);
VCFRecord record;

while (parser.next(record)) {
    for (const auto& sample : record.samples) {
        auto gt = sample.get_gt();
        auto dp = sample.get_dp();
        
        if (gt && dp) {
            std::cout << sample.name << ": " << *gt 
                      << " (depth=" << *dp << ")\n";
        }
    }
}
```

---

## Phase 10: Implementation Checklist

### 10.1 Core Features (MVP)

- [ ] Parse ##fileformat line
- [ ] Parse ##INFO lines
- [ ] Parse ##FORMAT lines
- [ ] Parse ##FILTER lines
- [ ] Parse #CHROM header line
- [ ] Parse data records (8 fixed fields)
- [ ] Handle missing values (".")
- [ ] Parse INFO field (key=value pairs)
- [ ] Parse FORMAT and sample columns
- [ ] Decode GT field
- [ ] Handle percent encoding/decoding
- [ ] Basic type validation
- [ ] Unit tests for all parsers

### 10.2 Extended Features
- [ ] Parse ##contig lines
- [ ] Parse ##ALT lines (symbolic alleles)
- [ ] Parse ##SAMPLE lines
- [ ] Parse ##PEDIGREE lines
- [ ] Validate Number=A,R,G correctly
- [ ] Handle variable-length genotype fields
- [ ] Support phased genotypes (PS field)
- [ ] Calculate genotype ordering for GL/PL
- [ ] Integration tests with real data

### 10.3 Advanced Features
- [ ] Structural variant support (##INFO SV keys)
- [ ] Breakend notation parsing
- [ ] gVCF support (<*> allele, END INFO)
- [ ] BCF binary format reading
- [ ] BCF binary format writing
- [ ] BGZF compression/decompression
- [ ] Tabix indexing support
- [ ] Random access by region
- [ ] Performance benchmarks

### 10.4 Quality Assurance
- [ ] Comprehensive error messages
- [ ] Input validation at all levels
- [ ] Memory leak testing
- [ ] Fuzz testing
- [ ] Thread safety analysis
- [ ] API documentation
- [ ] Usage examples
- [ ] Performance profiling

---

## Appendix A: Quick Reference

### A.1 Common Variant Types

| Type | REF | ALT | Example |
|------|-----|-----|---------|
| SNP | A | G | Single base change |
| Insertion | C | CAT | Insert AT after C |
| Deletion | CAT | C | Delete AT |
| MNP | AT | GC | Multi-nucleotide polymorphism |
| Complex | ATC | GG | Substitution + indel |

### A.2 Number Field Quick Reference

| Number | Meaning | Example (2 ALTs) |
|--------|---------|------------------|
| 0 | Flag | FLAG (no value) |
| 1 | Single value | DP=50 |
| 2 | Two values | HQ=10,20 |
| A | One per ALT | AC=5,3 (2 values) |
| R | One per allele (REF+ALTs) | AD=10,5,3 (3 values) |
| G | One per genotype | PL=0,10,100,15,20,90 (6 values for diploid) |
| . | Variable | ID=rs1;rs2;rs3 |

### A.3 GT Encoding Quick Reference

| GT | Meaning | Ploidy | Phased |
|----|---------|--------|--------|
| 0/0 | Homozygous REF | Diploid | No |
| 0/1 | Heterozygous | Diploid | No |
| 1/1 | Homozygous ALT1 | Diploid | No |
| 0\|1 | Heterozygous | Diploid | Yes |
| 0/1/2 | Three alleles | Triploid | No |
| 0 | REF | Haploid | N/A |
| 1 | ALT1 | Haploid | N/A |
| ./. | Missing | Diploid | No |

### A.4 Phred Scale Quick Reference

| Phred | Probability | Confidence |
|-------|-------------|------------|
| 10 | 10% error | 90% |
| 20 | 1% error | 99% |
| 30 | 0.1% error | 99.9% |
| 40 | 0.01% error | 99.99% |
| 50 | 0.001% error | 99.999% |
| 60 | 0.0001% error | 99.9999% |

**Formula:** `Phred = -10 * log10(P(error))`

---

## Appendix B: Error Handling Patterns

### B.1 Error Types

```cpp
enum class VCFError {
    // Parse errors
    INVALID_FORMAT,
    MISSING_FILEFORMAT,
    INVALID_HEADER,
    MALFORMED_LINE,
    INVALID_COLUMN_COUNT,
    
    // Type errors
    TYPE_MISMATCH,
    INVALID_INTEGER,
    INVALID_FLOAT,
    INVALID_GENOTYPE,
    
    // Validation errors
    UNDEFINED_INFO_KEY,
    UNDEFINED_FORMAT_KEY,
    UNDEFINED_FILTER,
    UNDEFINED_CONTIG,
    DUPLICATE_ID,
    
    // Data errors
    UNSORTED_POSITIONS,
    NON_CONTIGUOUS_CHROM,
    INVALID_REF_BASE,
    INVALID_ALT_BASE,
    
    // BCF specific
    INVALID_BCF_HEADER,
    UNSUPPORTED_BCF_VERSION,
    CORRUPTED_BGZF_BLOCK
};

class VCFException : public std::exception {
    VCFError error_;
    std::string message_;
    int line_number_;
public:
    VCFException(VCFError error, std::string message, int line = -1);
    const char* what() const noexcept override;
    VCFError error() const;
    int line() const;
};
```

### B.2 Error Reporting Best Practices

```cpp
// Good: Specific error with context
throw VCFException(
    VCFError::TYPE_MISMATCH,
    "INFO field 'DP' expects Integer, got 'abc'",
    line_number
);

// Good: Include position information
throw VCFException(
    VCFError::INVALID_ALT_BASE,
    "Invalid ALT base 'X' at chr1:12345",
    line_number
);

// Good: Suggest fix
throw VCFException(
    VCFError::UNDEFINED_INFO_KEY,
    "INFO key 'CUSTOM' not defined in header. Add ##INFO line.",
    line_number
);
```

---

## Appendix C: Performance Benchmarks

### C.1 Target Performance

**Typical File Sizes:**
- Small: 1,000 variants (~100 KB)
- Medium: 100,000 variants (~10 MB)
- Large: 10,000,000 variants (~1 GB)
- Very Large: 100,000,000+ variants (~10+ GB)

**Target Parse Speed:**
- Text VCF: 50-100 MB/s
- BCF: 200-500 MB/s
- With validation: 20-50% slower
- Streaming (no storage): Memory constant

---

## Final Notes

This specification covers VCF v4.3 comprehensively. Implement in phases:

1. **Phase 1**: Basic text parsing (days 1-3)
2. **Phase 2**: Complete validation (days 4-5)
3. **Phase 3**: Structural variants (day 6, optional)
4. **Phase 4**: gVCF support (day 7, optional)
5. **Phase 5**: BCF format (weeks 2-3, optional)

**Start with Phase 1, test thoroughly, then expand.**

Good luck with your implementation!</parameter>