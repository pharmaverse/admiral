echo "Waiting for rstudio-server to be ready.."

# Get rstudio url and open it
url=$(jq -r ".CODESPACE_NAME" /workspaces/.codespaces/shared/environment-variables.json)
url_rstudio="https://$url-8787.preview.app.github.dev"

# loop over url_rstudio until getting success status
MAX_RETRIES=10
WAIT_TIME=5

retries=0

while [ $retries -lt $MAX_RETRIES ]; do
  http_status=$(curl -s -o /dev/null -w "%{http_code}" "$url_rstudio")

  # Check status code (exit successfully if code in success/redirection category)
  if [ $http_status -ge 200 ]; then
    echo "Rstudio server ready ! You can now go to this url : $url_rstudio"
    exit 0
  else
    echo "$http_status"
    sleep $WAIT_TIME
    ((retries++))
  fi
done

echo "Error accessing rstudio ($url_rstudio) - Timeout reached."
exit 1