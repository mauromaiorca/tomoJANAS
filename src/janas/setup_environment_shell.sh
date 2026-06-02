#!/usr/bin/env bash
# ────────────────────────────────────────────────────
# janas shell integration helpers
# ────────────────────────────────────────────────────
# Usage: setup_environment_shell.sh [install|uninstall] [--shell SHELL] [--all] [--help]
#
#   install/uninstall    add or remove janas_activate_environment() / janas_deactivate_environment()
#   --shell SHELL        one of: bash, zsh, ksh, csh, tcsh, fish
#   --all                apply to every supported shell
#   -h, --help           show this message
#────────────────────────────────────────────────────

# File: setup_environment_shell.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology


# Determine project root (two levels up from this script)
JANAS_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# Virtualenv directory at project root
VENV="$JANAS_ROOT/.janas_env"

# POSIX snippet
posix_snippet='\
# >>> janas shell integration >>>
JANAS_ROOT="'"$JANAS_ROOT"'"
VENV="'"$VENV"'"
# janas_activate_environment: activate janas virtualenv
janas_activate_environment() {
  if [ -f "$VENV/bin/activate" ]; then
    . "$VENV/bin/activate"
  else
    echo "❌ janas virtualenv not found at $VENV"
  fi
}
# janas_deactivate_environment: deactivate any venv
deactivate() {
  if [ -n "$VIRTUAL_ENV" ]; then
    . "$VENV/bin/deactivate"
  else
    echo "❌ no active janas env"
  fi
}
# <<< janas shell integration <<<'

# fish snippet
fish_snippet='\
# >>> janas shell integration >>>
function janas_activate_environment
  if test -f "$VENV/bin/activate.fish"
    . "$VENV/bin/activate.fish"
  else
    echo "❌ janas virtualenv not found at $VENV"
  end
end
function janas_deactivate_environment
  if functions -q deactivate
    deactivate
  else
    echo "❌ no active janas env"
  end
end
# <<< janas shell integration <<<'

# Supported shells and their RC files
SHELLS=(bash zsh ksh csh tcsh fish)
RC_FILES=(
  "$HOME/.bashrc"
  "$HOME/.zshrc"
  "$HOME/.kshrc"
  "$HOME/.cshrc"
  "$HOME/.tcshrc"
  "$HOME/.config/fish/config.fish"
)

# Detect default shell
current_shell=$(basename "$SHELL" 2>/dev/null || echo bash)
default_shell=bash
for idx in "${!SHELLS[@]}"; do
  if [ "${SHELLS[idx]}" = "$current_shell" ]; then
    default_shell=$current_shell
    break
  fi
done

usage() {
  cat <<EOF
Usage: $0 [install|uninstall] [--shell SHELL] [--all] [-h|--help]

Commands:
  install       inject integration into rc file(s)
  uninstall     remove it

Options:
  --shell SHELL   target shell (choices: ${SHELLS[*]})
  --all           apply to all supported shells
  -h, --help      show this help and exit
EOF
}

# Parse action
if [[ "$1" == "install" || "$1" == "uninstall" ]]; then
  action=$1; shift
else
  action=install
fi

install_all=false
shell_choice=
# Parse flags
while (( "$#" )); do
  case "$1" in
    --shell) shell_choice=$2; shift 2 ;;
    --all)   install_all=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "❌ Unknown option: $1"; usage; exit 1 ;;
  esac
done

install_one() {
  rcf=$1; snippet=$2
  mkdir -p "$(dirname "$rcf")"
  if grep -F ">>> janas shell integration >>>" "$rcf" &>/dev/null; then
    echo "…already present in $rcf"
  else
    printf "\n%s\n" "$snippet" >>"$rcf"
    echo "✅ added snippet to $rcf"
  fi
}
uninstall_one() {
  rcf=$1
  if grep -F ">>> janas shell integration >>>" "$rcf" &>/dev/null; then
    sed -i.bak '/>>> janas shell integration >>>/,/<<< janas shell integration <<</d' "$rcf"
    echo "🗑️ removed snippet from $rcf (backup at $rcf.bak)"
  else
    echo "…nothing to remove in $rcf"
  fi
}

# Perform
case "$action" in
  install)
    if $install_all; then
      for idx in "${!SHELLS[@]}"; do
        snippet=$([[ "${SHELLS[idx]}" == "fish" ]] && echo "$fish_snippet" || echo "$posix_snippet")
        install_one "${RC_FILES[idx]}" "$snippet"
      done
    else
      [ -z "$shell_choice" ] && { shell_choice=$default_shell; echo "🔧 no --shell; using $default_shell"; }
      # find index
      idx=-1
      for i in "${!SHELLS[@]}"; do [ "${SHELLS[i]}" = "$shell_choice" ] && idx=$i && break; done
      if [ $idx -lt 0 ]; then echo "❌ unsupported shell: $shell_choice"; exit 1; fi
      snippet=$([[ "$shell_choice" == "fish" ]] && echo "$fish_snippet" || echo "$posix_snippet")
      install_one "${RC_FILES[idx]}" "$snippet"
    fi
    ;;
  uninstall)
    if $install_all; then
      for rcf in "${RC_FILES[@]}"; do uninstall_one "$rcf"; done
    else
      [ -z "$shell_choice" ] && { shell_choice=$default_shell; echo "🔧 no --shell; using $default_shell"; }
      idx=-1
      for i in "${!SHELLS[@]}"; do [ "${SHELLS[i]}" = "$shell_choice" ] && idx=$i && break; done
      if [ $idx -lt 0 ]; then echo "❌ unsupported shell: $shell_choice"; exit 1; fi
      uninstall_one "${RC_FILES[idx]}"
    fi
    ;;
  *) usage; exit 1 ;;
esac

