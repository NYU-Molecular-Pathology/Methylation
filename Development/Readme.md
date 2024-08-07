# Setup Methylation Developement Environment in RStudio

## Methylation Pipeline Development Directory
The methylation project can be opened in Rstudio in the following directory:
`/Volumes/CBioinformatics/Development/Methylation/pipeline/Methylation`
  + Make sure to **Git Pull** the latest version from Github before attempting to commit changes (under Tools > Version Control > Pull Branches)
  + If you do not have Rstudio configured with the Methylation repo the steps are outlined below

### 1. Configure git with Rstudio ############################################
  + If you do not have git installed, you can `brew install git` or download and install here: https://git-scm.com/download/
  + Install the `usethis` package if not already installed: `if(!require(usethis)){install.packages("usethis", dependencies=T)}`
  + Set your user name and email:
```R
usethis::use_git_config(user.name = "YourGithubName", user.email = "yourGithubEmail@mail.com")
```
### 1. Generate a Personal Access Token (PAT) for Authentication ############################################
To get a PAT from GitHub, simply call this function in RStudio:
```R
usethis::create_github_token()
```
  + This function will take you to github.com, assuming you’re already signed in.
  + You can also manage your personal access tokens from https://github.com/settings/tokens, by going to Settings ➡️ Developer settings ➡️ Personal access tokens.
  + You can click on “Generate new token” here, and adjust the Expiration behavior, including “No expiration”.<br/> 

<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/89248fc7b58dd5d0b3ea9a78ee88bf0219605b88/screenshots/new-personal-access-token-screenshot.png" width="80%"/><br/> 

  + Once you’re happy with the token’s Note, Expiration, and Scopes, click “Generate token”.  **You won’t be able to see this token again**,
so don’t close or navigate away from this browser window until you store the PAT locally.
  + Copy the PAT to the clipboard, anticipating what to do next: trigger a prompt that lets us store the PAT in the Git credential store.

If you use a password management app, such as 1Password or LastPass, you might want to also add this PAT (and its Note) to the entry for GitHub, 
where you’re already storing your username and password.

### 2. Set your personal access token in RStudio: ######################################################
If you don’t have gitcreds installed, install via `install.packages("gitcreds", dependencies=T)`. If you’ve already installed the latest usethis package, you will already have gitcreds, because the latest version of usethis includes gh and gh uses gitcreds.
  + Call `gitcreds::gitcreds_set()`. If you don’t have a PAT stored already, it will prompt you to enter your PAT.
  + **Paste your PAT from GitHub here.**
```
> gitcreds::gitcreds_set()

? Enter password or token: ghp_xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
-> Adding new credentials...
-> Removing credentials from cache...
-> Done.
```
Alternitavely:
```R
credentials::set_github_pat()
```
### 3. Restart RStudio and Verify settings ######################################################
```R
usethis::git_sitrep()
```
Your username and email should be stated correctly in the output. <br>
The report output shoud contain something like: `Personal access token: '<found in env var>'`
  + If you are still having troubles, read the output carefully. It might be that the PAT is not updated in your `.Renviron` file.
  + Call `usethis::edit_r_environ()` to update that file manually. You can reference additional information here: https://happygitwithr.com/https-pat.html

### 4. Create new RStudio Git Project ######################################################
  + In the RStudio Menubar, go to File ➡️ New Project... and select *Version Control* from the New Project Setup Wizard Menu
  + On the next page, select the option *Git* and paste the link to this repo in the *Repository URL* field on the following page: 
  https://github.com/NYU-Molecular-Pathology/Methylation/ </br>
  
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/2db0ce81ecb4bac1be64bf956f29ba16f430bc9c/screenshots/new-git-project-rstudio.png" width="50%"/><br/> 

  + It should automatically create a directory for you named Methylation. Check the box "Open in new session" and click "Create Project"
  + This will clone the latest Methylation repository files locally in a new RStudio project environment
  + If you want to load the existing Development project, in Rstudio open a Project in a new session from `/Volumes/CBioinformatics/Development/Methylation/Pipeline/Methylation` and then be sure to Git > Pull the latest updates first
