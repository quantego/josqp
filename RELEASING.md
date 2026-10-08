# Releasing jOSQP

Releases go to Maven Central through the Sonatype Central Portal. Publishing is
fully automated by GitHub Actions: pushing a `v*` tag runs
`.github/workflows/release.yml`, which builds and tests the artifact, signs it
and uploads it. `autoPublish=true` in `pom.xml` means the deployment is
published right away, so a green run ends with the version live in Central.

**Secrets never belong in the repository, in a ticket or in a chat.** They live
in the GitHub repository secrets and, if you also want to release locally, in
the GPG keyring and `~/.m2/settings.xml` of your machine.

## One-time setup

### 1. GPG signing key

Central rejects unsigned artifacts. Generate a key once:

```
gpg --quick-generate-key "jOSQP Release <mail@example.com>" rsa4096 sign never
KEYID=$(gpg --list-secret-keys --with-colons | awk -F: '/^sec:/ { print $5; exit }')
gpg --keyserver keyserver.ubuntu.com --send-keys "$KEYID"
```

`--quick-generate-key` asks for the passphrase in the normal pinentry dialog and
never puts it on the command line. The wizard variant `gpg --full-generate-key`
does the same thing with more questions. Instead of a keyserver the public key
can also be uploaded in the Central Portal account.

### 2. Central Portal token

Account at <https://central.sonatype.com>, namespace `com.quantego` (already
verified: `com.quantego:josqp:0.6.5` was published from it). Under *Account ->
Generate User Token* you get a username and a password.

### 3. GitHub repository secrets

| Secret | Value |
|---|---|
| `CENTRAL_USERNAME` | user-token username |
| `CENTRAL_PASSWORD` | user-token password |
| `GPG_PRIVATE_KEY` | ASCII-armoured private key from step 1 |
| `GPG_PASSPHRASE` | passphrase of that key |

The simplest route needs no GitHub token at all: create the four secrets by hand
under <https://github.com/quantego/josqp/settings/secrets/actions/new>. The
private key can be copied straight to the clipboard with
`gpg --armor --export-secret-keys "$KEYID" | pbcopy`.

With the GitHub CLI (`gh auth login --scopes repo` once, `KEYID` from step 1):

```
R=quantego/josqp
gpg --armor --export-secret-keys "$KEYID" | gh secret set GPG_PRIVATE_KEY --repo "$R"
printf "GPG-Passphrase: " && read -rs P && echo && printf '%s' "$P" | gh secret set GPG_PASSPHRASE --repo "$R"
printf "Central-Token Benutzername: " && read -r CU && printf '%s' "$CU" | gh secret set CENTRAL_USERNAME --repo "$R"
printf "Central-Token Passwort: " && read -rs CP && echo && printf '%s' "$CP" | gh secret set CENTRAL_PASSWORD --repo "$R"
gh secret list --repo "$R"
```

Nothing is echoed, stored in a file or written to the shell history. Two zsh
traps to avoid when pasting these lines: an inline `# comment` after a command
is read as an *argument* unless `interactive_comments` is set, and an unquoted
`var?prompt` argument to `read` is treated as a glob pattern.

Instead of the browser flow a **fine-grained personal access token** works too. It
is created in your **personal** account (*Settings -> Developer settings ->
Personal access tokens -> Fine-grained tokens*), but `quantego` has to be picked
as the **resource owner**, because the repository belongs to the organization;
the organization must allow that first (*Organization -> Settings -> Personal
access tokens -> Allow access via fine-grained personal access tokens*), and if
approval is required an organization owner has to approve the token. Then set
*Repository access* to "Only select repositories" and pick `quantego/josqp`.
Repository permissions:

| Permission | Why |
|---|---|
| `Secrets: Read and write` | set the four repository secrets |
| `Metadata: Read` | mandatory, added automatically |
| `Contents: Read and write` | only if the same token should push branch and tag |
| `Workflows: Read and write` | only for pushes that touch `.github/workflows/` |

`gh`'s own help recommends the `GH_TOKEN` environment variable for fine-grained
tokens rather than `--with-token`:

```
printf "GitHub-Token: " && read -rs T && echo && export GH_TOKEN="$T"
gh auth status && gh secret list --repo quantego/josqp
```

## Release

1. Set `<version>` in `pom.xml`, move the `## [Unreleased]` entries of
   [CHANGELOG.md](CHANGELOG.md) into a dated section (the date is the day you publish) and commit.
2. Wait for portable CI on Java 8, 17 and 21, and the Vector CI on Java 26. Confirm the
   release tag equals the POM version; the workflow rejects mismatches and manual branch runs.
   Only then tag and push after the explicit publication approval:

   ```
   git tag -a v1.0.0 -m "jOSQP 1.0.0"
   git push origin master v1.0.0
   ```

3. Watch the *Release* workflow. Afterwards create the GitHub release from the
   tag and use the matching changelog section as its notes.
4. Open a fresh `## [Unreleased]` section for the next cycle.

Do **not** add `-Pvector` to the release: the Vector API backend needs JDK 26
and the `jdk.incubator.vector` module and must not end up in the plain Java 8
artifact.

## Releasing locally instead

`mvn -Prelease deploy` does the same on a machine that has the keyring and the
Central token; the token then belongs in `~/.m2/settings.xml` (permissions `600`):

```xml
<settings>
  <servers>
    <server>
      <id>central</id>
      <username>TOKEN-USERNAME</username>
      <password>TOKEN-PASSWORD</password>
    </server>
  </servers>
</settings>
```

The passphrase comes from `MAVEN_GPG_PASSPHRASE` or a running `gpg-agent`.

## Verifying a release locally

```bash
mvn test                            # portable suite
mvn -Pvector test                   # optional suite, JDK 26
mvn -Prelease -Dgpg.skip=true verify # tests + all three jars, without a signing key
```

The local unsigned verification does not test GPG signing or Central authentication. Those
require the configured repository secrets. Never enable the Vector profile for a Central release.
