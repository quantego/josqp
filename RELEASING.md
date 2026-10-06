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

Central rejects unsigned artifacts. Generate a key once (interactive, the
passphrase is asked for with hidden input and never ends up in the shell
history):

```
gpg --full-generate-key            # RSA 4096, no expiry is fine
gpg --list-secret-keys --keyid-format=long
gpg --keyserver keyserver.ubuntu.com --send-keys <KEY-ID>
```

Then copy the private key into the repository secret `GPG_PRIVATE_KEY`:

```
gpg --armor --export-secret-keys <KEY-ID> | pbcopy
# GitHub -> Settings -> Secrets and variables -> Actions -> New repository secret
```

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

## Release

1. Set `<version>` in `pom.xml`, move the `## [Unreleased]` entries of
   [CHANGELOG.md](CHANGELOG.md) into a dated section and commit.
2. Tag and push:

   ```
   git tag -a v1.0.0 -m "jOSQP 1.0.0"
   git push origin master --tags
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
mvn test                            # 90 tests
mvn -Pvector test                   # 93 tests, JDK 26
mvn -Prelease -DskipTests package   # main, sources and javadoc jars in target/
```
