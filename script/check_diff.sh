#!/usr/bin/env bash

# Source and destination paths
srcPath=$1
dstPath=$2
commandPath=$3
jediPrefix=$4

# Initial check
if test -f "${dstPath}"; then
  # Destination file exists
  if cmp -s ${srcPath} ${dstPath}; then
    # Destination file and source file are similar, exiting
    exit 0
  fi
fi

# Prepare destination directory
dstDir=$(dirname "${dstPath}")
mkdir -p ${dstDir}

# Copy file
cp -f ${srcPath} ${dstPath}.tmp

# Get file extension
srcFile=$(basename "${srcPath}")
srcExt="${srcFile##*.}"


if test "${srcExt}" = "h" -o "${srcExt}" = "cc"; then
  if test "${jediPrefix}" = "ON"; then
    # Update content
    if test "${srcFile}" = "Variables.h"; then
      sed -i -e s/"Variables"/"JediVariables"/g ${dstPath}.tmp
      sed -i -e s/"JediVariablesBase"/"VariablesBase"/g ${dstPath}.tmp
    elif test "${srcFile}" = "Variables.cc"; then
      sed -i -e s/"Variables"/"JediVariables"/g ${dstPath}.tmp
      sed -i -e s/"JediVariablesBase"/"VariablesBase"/g ${dstPath}.tmp
      sed -i -e s/"\/JediVariables"/"\/Variables"/g ${dstPath}.tmp
    else
      sed -i -e s/"Variables::Variables"/"JediVariables::JediVariables"/g ${dstPath}.tmp
      sed -i -e s/"oops::Variables"/"oops::JediVariables"/g ${dstPath}.tmp
      sed -i -e s/" Variables"/" JediVariables"/g ${dstPath}.tmp
      sed -i -e s/"<Variables>"/"<JediVariables>"/g ${dstPath}.tmp
      sed -i -e s/"^Variables"/"JediVariables"/g ${dstPath}.tmp
    fi
  elif test "${jediPrefix}" = "OFF"; then
    sed -i -e s/"oops::Variables"/"Variables"/g ${dstPath}.tmp
  else
    echo "Wrong JEDIPREFIX: "${JEDIPREFIX}
    exit 1
  fi
fi

# Check patch size
if test -f "${dstPath}.patch"; then
  patchLength=`cat ${dstPath}.patch | wc -l`
  if test "${patchLength}" = "0"; then
    rm -f ${dstPath}.patch
  fi
fi

if test -f "${dstPath}.patch"; then
  # Apply residual patch
  cp -f ${dstPath}.tmp ${dstPath}.tmp.bak
  patch -s ${dstPath}.tmp ${dstPath}.patch

  # Compare and update if needed
  if cmp -s ${dstPath}.tmp ${dstPath}; then
    rm -f ${dstPath}.tmp ${dstPath}.tmp.bak
  else
    echo "--  - Update needed for: ${dstPath}"
    echo "echo \"Difference between ${dstPath}.tmp and ${dstPath}: [u]pdate, [m]eld or keep (any key)?\"; read action; if test \"\${action}\" = \"u\"; then cp -f ${dstPath}.tmp ${dstPath}; elif test \"\${action}\" = \"m\"; then meld ${dstPath}.tmp ${dstPath} ${dstPath}.tmp.bak; fi; diff -u ${dstPath}.tmp.bak ${dstPath} > ${dstPath}.patch; rm -f ${dstPath}.tmp ${dstPath}.tmp.bak" >> ${commandPath}
    echo "" >> ${commandPath}
 fi
else
  if test -f "${dstPath}"; then
    # Create patch if needed
    diff -u ${dstPath}.tmp ${dstPath} > ${dstPath}.patch
    patchLength=`cat ${dstPath}.patch | wc -l`
    if test "${patchLength}" = "0"; then 
      rm -f ${dstPath}.patch
      rm -f ${dstPath}.tmp
    else
      echo "--  - New patch needed for: "${dstPath}
      rm -f ${dstPath}.patch
      echo "echo \"Difference between ${dstPath}.tmp and ${dstPath}: [u]pdate, [m]eld or keep (any key)?\"; read action; if test \"\${action}\" = \"u\"; then cp -f ${dstPath}.tmp ${dstPath}; elif test \"\${action}\" = \"m\"; then meld ${dstPath}.tmp ${dstPath}; fi; diff -u ${dstPath}.tmp ${dstPath} > ${dstPath}.patch; patchLength=\`cat ${dstPath}.patch | wc -l\`; if test "\$\{patchLength\}" = "0"; then rm -f ${dstPath}.patch;fi;rm -f ${dstPath}.tmp" >> ${commandPath}
      echo "" >> ${commandPath}
    fi
  else
    # New file
    mv ${dstPath}.tmp ${dstPath}
  fi
fi
