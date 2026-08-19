#!/bin/bash

# Fortran ファイル整形スクリプト (変更があったファイルのみ書き換え)
# 使用法: format_fortran_script_update_only.sh

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

# 処理したファイル数をカウント
processed=0
errors=0

shopt -s nullglob  # マッチしない場合は空の配列
for file in *.f90; do
    if [ -f "$file" ]; then
        # findent で整形 → 一時ファイルに保存
        if findent -i2 -Rr < "$file" > "${file}.tmp"; then
            # 差分を比較
            if ! cmp -s "$file" "${file}.tmp"; then
                # バックアップ作成
                backup_file="$BACKUP_DIR/$(basename "$file")"
                cp "$file" "$backup_file"

                mv "${file}.tmp" "$file"
                echo "更新: $file"
                processed=$((processed + 1))
            else
                rm -f "${file}.tmp"
            fi
        else
            echo "  ✗ エラー: findent の実行に失敗しました ($file)"
            rm -f "${file}.tmp"
            errors=$((errors + 1))
        fi
    fi
done

echo "完了: 更新=$processed, エラー=$errors"
