#!/usr/bin/env python
"""
rna_structure_parser.py
-----------------------
RNA構造データの解析ユーティリティ

サポートされている形式:
- dbn形式 (dot-bracket notation)
- bpseq形式 (base pair sequence)

"""
from __future__ import annotations
from typing import List, Tuple, Dict, Optional, Union
import numpy as np
from pathlib import Path


class RNAStructureParser:
    """RNA構造データの汎用パーサークラス"""
    
    def __init__(self, filepath: Optional[str] = None, sequence: str = None, structure: str = None, name: str = None):
        """
        初期化
        
        Args:
            filepath: ファイルパス（dbnまたはbpseq）
            sequence: RNA配列（直接指定する場合）
            structure: 構造データ（直接指定する場合）
            name: RNA名（直接指定する場合）
        """
        self.filepath = filepath
        self.name = name
        self.length = None
        self.page_number = None
        self.sequence = sequence
        self.structure = structure
        self.pairs = None
        
        if filepath:
            if filepath.endswith('.dbn'):
                self._parse_dbn_file()
            elif filepath.endswith('.bpseq'):
                self._parse_bpseq_file()
            else:
                raise ValueError(f"Unsupported file format: {filepath}")
        elif sequence and structure:
            self.length = len(sequence)
            self.pairs = self._extract_pairs_from_dotbracket(structure)
        else:
            raise ValueError("Either filepath or sequence+structure must be provided")
    
    def _parse_dbn_file(self) -> None:
        """dbnファイルを解析してデータを抽出"""
        with open(self.filepath, 'r') as f:
            lines = [line.strip() for line in f.readlines()]
        
        # 空行を除去
        lines = [line for line in lines if line]
        
        if len(lines) < 5:
            raise ValueError(f"Invalid dbn file format: {self.filepath}")
        
        # メタデータを解析
        self.name = lines[0].split(':')[1].strip() if ':' in lines[0] else None
        self.length = int(lines[1].split(':')[1].strip()) if ':' in lines[1] else None
        self.page_number = int(lines[2].split(':')[1].strip()) if ':' in lines[2] else None
        
        # 配列と構造を取得
        self.sequence = lines[3]
        self.structure = lines[4]
        
        # 長さの整合性をチェック
        if len(self.sequence) != len(self.structure):
            raise ValueError(f"Sequence and structure length mismatch in {self.filepath}")
        
        if self.length and len(self.sequence) != self.length:
            raise ValueError(f"Declared length {self.length} doesn't match actual length {len(self.sequence)} in {self.filepath}")
        
        # 塩基対を抽出
        self.pairs = self._extract_pairs_from_dotbracket(self.structure)
    
    def _parse_bpseq_file(self) -> None:
        """bpseqファイルを解析してデータを抽出"""
        sequence = []
        pairs = []
        
        with open(self.filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split()
                if len(parts) >= 3:
                    pos = int(parts[0])
                    base = parts[1]
                    partner = int(parts[2]) if parts[2] != '0' else 0
                    
                    sequence.append(base)
                    if partner > pos:  # i < j のペアのみ記録
                        pairs.append((pos, partner))
        
        self.sequence = ''.join(sequence)
        self.length = len(self.sequence)
        self.pairs = pairs
        
        # dot-bracket記法も生成
        self.structure = self._generate_dotbracket()
        
        # ファイル名からRNA名を推定
        if not self.name:
            self.name = Path(self.filepath).stem
    
    def _extract_pairs_from_dotbracket(self, structure: str) -> List[Tuple[int, int]]:
        """dot-bracket記法から塩基対リストを抽出（1-indexedで返す）"""
        stack = []
        pairs = []
        
        for i, char in enumerate(structure):
            if char == '(':
                stack.append(i + 1)  # 1-indexedに変換
            elif char == ')':
                if stack:
                    j = stack.pop()
                    pairs.append((j, i + 1))  # (開く位置, 閉じる位置)
        
        # i < j の順序でソート
        pairs.sort()
        return pairs
    
    def _generate_dotbracket(self) -> str:
        """塩基対リストからdot-bracket記法を生成"""
        structure = ['.'] * self.length
        
        for i, j in self.pairs:
            structure[i-1] = '('  # 0-indexedに変換
            structure[j-1] = ')'
        
        return ''.join(structure)
    
    def to_bpseq_format(self) -> List[Tuple[int, str, int]]:
        """bpseq形式のデータを生成（位置、塩基、ペア相手）"""
        # ペア情報を辞書に変換
        pair_dict = {}
        for i, j in self.pairs:
            pair_dict[i] = j
            pair_dict[j] = i
        
        bpseq_data = []
        for i, base in enumerate(self.sequence, 1):
            partner = pair_dict.get(i, 0)
            bpseq_data.append((i, base, partner))
        
        return bpseq_data
    
    def get_structure_info(self) -> Dict:
        """構造情報を辞書形式で返す"""
        return {
            "name": self.name,
            "length": len(self.sequence) if self.sequence else 0,
            "page_number": self.page_number,
            "sequence": self.sequence,
            "structure": self.structure,
            "num_pairs": len(self.pairs) if self.pairs else 0,
            "pairs": self.pairs
        }
    
    def save_as_fasta(self, output_path: str) -> None:
        """FASTA形式で保存"""
        with open(output_path, 'w') as f:
            f.write(f">{self.name}\n")
            f.write(f"{self.sequence}\n")
    
    def save_as_bpseq(self, output_path: str) -> None:
        """bpseq形式で保存"""
        bpseq_data = self.to_bpseq_format()
        with open(output_path, 'w') as f:
            for pos, base, partner in bpseq_data:
                f.write(f"{pos} {base} {partner}\n")
    
    def save_as_dbn(self, output_path: str) -> None:
        """dbn形式で保存"""
        with open(output_path, 'w') as f:
            f.write(f"#Name: {self.name}\n")
            f.write(f"#Length: {self.length}\n")
            if self.page_number:
                f.write(f"#PageNumber: {self.page_number}\n")
            f.write(f"{self.sequence}\n")
            f.write(f"{self.structure}\n")


def parse_rna_files(data_dir: str, file_pattern: str = "*.dbn", max_files: Optional[int] = None) -> List[RNAStructureParser]:
    """
    指定ディレクトリ内のRNAファイルを一括解析
    
    Args:
        data_dir: ファイルが格納されているディレクトリ
        file_pattern: ファイルパターン（*.dbn, *.bpseq など）
        max_files: 処理する最大ファイル数（テスト用）
    
    Returns:
        RNAStructureParserオブジェクトのリスト
    """
    data_dir = Path(data_dir)
    files = list(data_dir.glob(file_pattern))
    
    if max_files:
        files = files[:max_files]
    
    parsers = []
    for file_path in files:
        try:
            parser = RNAStructureParser(str(file_path))
            parsers.append(parser)
        except Exception as e:
            print(f"Error parsing {file_path}: {e}")
    
    return parsers


# 後方互換性のためのエイリアス
DbnParser = RNAStructureParser


if __name__ == "__main__":
    # テスト用
    dbn_dir = "/Users/ootagakitakumi/Library/Mobile Documents/com~apple~CloudDocs/大学院/浅井研究室/研究/lincapr/data/processed/bprna/dbn_pagenumber1"
    
    # 1つのファイルをテスト
    test_file = f"{dbn_dir}/bpRNA_CRW_1236.dbn"
    parser = RNAStructureParser(test_file)
    
    print("=== Test RNA structure parsing ===")
    print(f"Name: {parser.name}")
    print(f"Length: {parser.length}")
    print(f"Actual length: {len(parser.sequence)}")
    print(f"Number of pairs: {len(parser.pairs)}")
    print(f"First 5 pairs: {parser.pairs[:5]}")
    
    # FASTA形式で保存テスト
    parser.save_as_fasta("/tmp/test.fa")
    parser.save_as_bpseq("/tmp/test.bpseq")
    print("Test files saved to /tmp/test.fa and /tmp/test.bpseq")
