# Setup Methylation Developement Environment in RStudio
### 1. Install git ##########################################################
If you do not have git installed, you can download and install here: https://git-scm.com/download/

### 2. Configure git with Rstudio ############################################

#### Install usethis package if not already installed:
```R
if(!require(usethis)){install.packages("usethis", dependencies=T)}
```
#### Set your user name and email:
```ruby
usethis::use_git_config(user.name = "YourGithubName", user.email = "yourGithubEmail@mail.com")
```
#### Generate a Personal Access Token (PAT) for Authentication:
To get a PAT from GitHub, simply call this function from RStudio:
```R
usethis::create_github_token()
```
  + This function will take you to github.com, assuming you’re already signed in, you can manage your personal access tokens from https://github.com/settings/tokens, 
also by going to Settings ➡️ Developer settings ➡️ Personal access tokens. You can click on “Generate new token” here, and adjust the Expiration behaviour as you see fit, including “No expiration”. 
  + Once you’re happy with the token’s Note, Expiration, and Scopes, click “Generate token”.  **You won’t be able to see this token again**,
so don’t close or navigate away from this browser window until you store the PAT locally.
  + Copy the PAT to the clipboard, anticipating what to do next: trigger a prompt that lets us store the PAT in the Git credential store.

If you use a password management app, such as 1Password or LastPass, you might want to also add this PAT (and its Note) to the entry for GitHub, 
where you’re already storing your username and password.

#### Set your personal access token in RStudio:
If you don’t have gitcreds installed, install via install.packages("gitcreds"). <br>
If you’ve already installed the latest usethis package, you will already have gitcreds, because the newest usethis uses gh and gh uses gitcreds.

Call gitcreds::gitcreds_set(). If you don’t have a PAT stored already, it will prompt you to enter your PAT. Paste your PAT from GitHub here.
```bash
> gitcreds::gitcreds_set()

? Enter password or token: ghp_xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
-> Adding new credentials...
-> Removing credentials from cache...
-> Done.
```
Alternitavely:
```R
credentials::set_github_pat("YourPAT")
```
#### 4. Restart RStudio ######################################################
#### 5. Verify settings ######################################################

```R
usethis::git_sitrep()
```
Your username and email should be stated correctly in the output. <br>
The report output shoud contain something like:
`Personal access token: '<found in env var>'`

  + If you are still having troubles, read the output carefully.
  + It might be that the PAT is not updated in your `.Renviron` file.
  + Call `usethis::edit_r_environ()` to update that file manually.
  + You can reference additional information here: https://happygitwithr.com/https-pat.html
