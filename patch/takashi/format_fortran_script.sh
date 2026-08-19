#!/bin/bash

# Fortran ファイル整形スクリプト (確実版)
# 使用法: ./format_fortran.sh

set -e

TARGET_DIR="."

# findent の存在確認
if ! command -v findent &> /dev/null; then
    echo "エラー: findent が見つかりません"
    echo "インストール方法:"
    echo "  Ubuntu/Debian: sudo apt-get install findent"
    echo "  CentOS/RHEL:   sudo yum install findent"
    echo "  macOS:         brew install findent"
    exit 1
fi

# バックアップディレクトリの作成
BACKUP_DIR="${TARGET_DIR}/backup_$(date +%Y%m%d_%H%M%S)"
echo "バックアップディレクトリを作成: $BACKUP_DIR"
mkdir -p "$BACKUP_DIR"

echo "対象ディレクトリ: $TARGET_DIR"
echo "Fortran ファイルを検索中..."

# 処理したファイル数をカウント
processed=0
errors=0

# 直接的なアプローチ: forループでglob展開を使用
shopt -s nullglob  # マッチしない場合は空の配列
for file in *.f90; do
    if [ -f "$file" ]; then
        echo "処理中: $file"
        
        # バックアップ作成
        backup_file="$BACKUP_DIR/$(basename "$file")"
        cp "$file" "$backup_file"
        
        # findent で整形
        if findent -i2 -Rr < "$file" > "${file}.tmp"; then
            mv "${file}.tmp" "$file"
            echo "  ✓ 完了"
            processed=$((processed + 1))
        else
            echo "  ✗ エラー: findent の実行に失敗しました"
            rm -f "${file}.tmp"
            errors=$((errors + 1))
        fi
    fi
done

echo ""
echo "処理完了!"
echo "処理済みファイル: $processed"
echo "エラー: $errors"
echo "バックアップ場所: $BACKUP_DIR"

if [ "$errors" -gt 0 ]; then
    echo ""
    echo "警告: エラーが発生したファイルがあります。"
    echo "必要に応じてバックアップから復元してください。"
    exit 1
fi