# Releasing jOSQP

Maintainer checklist for publishing `com.quantego:josqp` to Maven Central.
The [release workflow](.github/workflows/release.yml) builds, tests, signs and deploys
the portable Java 8 artifact when a `v*` tag is pushed. `autoPublish=true` means
a successful deployment is published automatically.

## Credentials

Configure these [repository Actions secrets](https://github.com/quantego/josqp/settings/secrets/actions):

| Secret | Value |
|---|---|
| `CENTRAL_USERNAME` | Central Portal user-token username |
| `CENTRAL_PASSWORD` | Central Portal user-token password |
| `GPG_PRIVATE_KEY` | ASCII-armoured private signing key |
| `GPG_PASSPHRASE` | Passphrase for that key |

Publish the corresponding public key on a Central-supported keyserver, such as
`keys.openpgp.org`; complete its email verification and confirm the key is retrievable.
No public-key secret is needed. Keep the private key, passphrase and revocation
certificate securely backed up. Never commit credentials or paste them into issues.

See [Sonatype's PGP requirements](https://central.sonatype.org/publish/requirements/gpg/)
and [keys.openpgp.org's upload instructions](https://keys.openpgp.org/about/usage-gnupg/)
for key setup, renewal and publication.

## Release checklist

1. Set the new version in `pom.xml`, update installation examples and move the
   relevant [changelog](CHANGELOG.md) entries into a dated release section.
2. Verify the portable release artifacts locally, then commit and push the branch:

   ```bash
   mvn -Prelease -Dgpg.skip=true verify
   ```

   This checks tests and packaging, but does not test signing or Central authentication.
   Check that the jars contain the license/notice texts and only Java 8 library classes.
3. Wait for [CI](.github/workflows/ci.yml) on Java 8, 17 and 21 and the optional
   Vector suite on Java 26 to pass for that commit.
4. Create and push an annotated tag `vX.Y.Z` matching the POM version.
   **Pushing the tag starts automatic publication.** Do not move an existing published tag.
5. Confirm the [Release job](https://github.com/quantego/josqp/actions/workflows/release.yml)
   succeeded and the version is available on Maven Central. Create the GitHub release
   with the matching changelog notes and open a new `Unreleased` section.

The release runs on JDK 21 without `-Pvector`; the experimental Vector backend is
not part of the published portable jar.
