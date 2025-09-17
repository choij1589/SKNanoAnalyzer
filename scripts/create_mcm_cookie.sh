#!/bin/bash
# Helper script to create MCM cookie file from browser cookie

echo "Please paste the cookie value from your browser (Cookie header value):"
read -r cookie_value

# Create the cookie file in a simple format
mkdir -p ~/private
echo "$cookie_value" > ~/private/prod-cookie.txt
chmod 600 ~/private/prod-cookie.txt

echo "Cookie saved to ~/private/prod-cookie.txt"
echo "You can now use the fetchMCMFilterEfficiency.py script"
